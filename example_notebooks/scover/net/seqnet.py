import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import pytorch_lightning as pl
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks import TQDMProgressBar
from ray import tune
from ray.tune.integration.pytorch_lightning import TuneReportCallback
from ray.tune.search.hebo import HEBOSearch
from ray.tune import CLIReporter, JupyterNotebookReporter
from ray.tune.schedulers import ASHAScheduler, PopulationBasedTraining


def loguniform(a=0, b=1, size=None):
    r"""
    Returns random samples from loguniform distribution 
      between 'a' and 'b' of a particular size
    """
    X = 10**((np.log10(b) - np.log10(a)) * np.random.random_sample(size) +
             np.log10(a))
    return X


def train_scover_bs(config,
                    train_data,
                    val_data,
                    max_epochs=30,
                    num_gpus=1,
                    use_elu=False,
                    use_smooth_l1=False):
    if num_gpus > 0:
        assert torch.cuda.is_available(
        ), "num_gpus > 0 but CUDA is not available"
    train_loader = DataLoader(train_data,
                              batch_size=config["batch_size"],
                              num_workers=0,
                              shuffle=True)
    val_loader = DataLoader(val_data,
                            batch_size=config["batch_size"],
                            num_workers=0,
                            shuffle=False)
    seq_length = train_loader.dataset[0][0].shape[1]
    output_size = train_loader.dataset[0][1].shape[0]
    sn = SeqNet(seq_length=seq_length,
                output_size=output_size,
                learning_rate=config["learning_rate"],
                motif_length=15,
                num_motifs=600,
                sigma_motifs=config["sigma_motifs"],
                sigma_net=config["sigma_net"],
                use_elu=use_elu,
                use_smooth_l1=use_smooth_l1)
    tune_callback = TuneReportCallback({
        "loss": "ptl/val_loss",
    },
                                       on="validation_end")
    pb_callback = TQDMProgressBar(refresh_rate=0)
    trainer = pl.Trainer(max_epochs=max_epochs,
                         gpus=num_gpus,
                         callbacks=[tune_callback, pb_callback],
                         logger=TensorBoardLogger(
                             save_dir=tune.get_trial_dir(),
                             name="",
                             version="."))
    trainer.fit(sn, train_dataloaders=train_loader, val_dataloaders=val_loader)


def train_scover(config,
                 train_loader,
                 val_loader,
                 max_epochs=30,
                 num_gpus=1,
                 use_elu=False,
                 use_smooth_l1=False):
    if num_gpus > 0:
        assert torch.cuda.is_available(
        ), "num_gpus > 0 but CUDA is not available"
    seq_length = train_loader.dataset[0][0].shape[1]
    output_size = train_loader.dataset[0][1].shape[0]
    sn = SeqNet(seq_length=seq_length,
                output_size=output_size,
                learning_rate=config["learning_rate"],
                motif_length=15,
                num_motifs=600,
                sigma_motifs=config["sigma_motifs"],
                sigma_net=config["sigma_net"],
                use_elu=use_elu,
                use_smooth_l1=use_smooth_l1)
    tune_callback = TuneReportCallback({
        "loss": "ptl/val_loss",
    },
                                       on="validation_end")
    trainer = pl.Trainer(max_epochs=max_epochs,
                         gpus=num_gpus,
                         progress_bar_refresh_rate=0,
                         callbacks=[tune_callback],
                         logger=TensorBoardLogger(
                             save_dir=tune.get_trial_dir(),
                             name="",
                             version="."))
    trainer.fit(sn, train_dataloaders=train_loader, val_dataloaders=val_loader)


def tune_scover_asha_hyperopt(train_data,
                              val_data,
                              num_samples=10,
                              num_epochs=50,
                              gpus_per_trial=1,
                              local_dir="~/ray_results",
                              prefix=None,
                              use_elu=False,
                              use_smooth_l1=False):
    r"""
    Experimental method for hyperparameter search using 
      hyperopt + Ray Tune
    """
    import numpy as np
    from hyperopt import hp
    from ray.tune.search.hyperopt import HyperOptSearch
    if prefix != None:
        run_name = prefix + "_tune_scover"
    else:
        run_name = "tune_scover"
    space = {
        "learning_rate": hp.loguniform("learning_rate", np.log(1e-5),
                                       np.log(1e-2)),
        "sigma_motifs": hp.loguniform("sigma_motifs", np.log(1e-7),
                                      np.log(1e-3)),
        "sigma_net": hp.loguniform("sigma_net", np.log(1e-5), np.log(1e-2)),
        "batch_size": hp.choice("batch_size", [64, 128, 256, 512])
    }
    scheduler = ASHAScheduler(metric='loss',
                              mode='min',
                              max_t=num_epochs,
                              grace_period=1,
                              reduction_factor=2)
    hyperopt_search = HyperOptSearch(space, metric='loss', mode='min')
    reporter = JupyterNotebookReporter(
        overwrite=True,
        parameter_columns=[
            "batch_size", "learning_rate", "sigma_motifs", "sigma_net"
        ],
        metric_columns=["loss", "training_iteration"])
    analysis = tune.run(tune.with_parameters(train_scover_bs,
                                             train_data=train_data,
                                             val_data=val_data,
                                             max_epochs=num_epochs,
                                             num_gpus=gpus_per_trial,
                                             use_elu=use_elu,
                                             use_smooth_l1=use_smooth_l1),
                        resources_per_trial={
                            "cpu": 1,
                            "gpu": gpus_per_trial
                        },
                        search_alg=hyperopt_search,
                        num_samples=num_samples,
                        scheduler=scheduler,
                        progress_reporter=reporter,
                        name=run_name,
                        max_failures=3,
                        local_dir=local_dir)
    print("Best hyperparameters found were: ",
          analysis.get_best_config(metric='loss', mode='min'))
    return analysis


class SeqNet(pl.LightningModule):
    r"""
    EXPERIMENTAL
    """

    def __init__(self,
                 seq_length,
                 output_size,
                 sigma_motifs=1e-5,
                 sigma_net=1e-3,
                 num_motifs=600,
                 motif_length=12,
                 beta_1_min=1e-15,
                 beta_1_max=1e-3,
                 beta_2_min=1e-10,
                 beta_2_max=1e-3,
                 learning_rate=5e-4,
                 use_elu=False,
                 use_smooth_l1=False):
        super(SeqNet, self).__init__()
        # Convolution and pooling
        self.conv_1 = nn.Conv2d(1, num_motifs, (motif_length, 4), bias=False)
        self.pool_1 = nn.MaxPool2d(((seq_length - motif_length + 1), 1))
        self.nonlin = nn.ELU() if use_elu else nn.ReLU()
        # Fully-connected layer
        self.fc = nn.Linear(num_motifs, output_size, bias=True)
        self.sigma_motifs = sigma_motifs
        self.sigma_net = sigma_net
        self.num_motifs = num_motifs
        self.motif_length = motif_length
        self.seq_length = seq_length
        self.output_size = output_size
        self.learning_rate = learning_rate
        self.beta = [
            loguniform(beta_1_min, beta_1_max),
            loguniform(beta_2_min, beta_2_max)
        ]
        self.apply(self.weights_init)
        if use_smooth_l1:
            self.loss = nn.SmoothL1Loss()
        else:
            self.loss = nn.MSELoss()

    def weights_init(self, m):
        r"""
        Method for initializing weights given `sigma_motifs`
          and `sigma_net`
        """
        classname = m.__class__.__name__
        if classname.find('Conv2d') != -1:
            m.weight.data.normal_(0.0, self.sigma_motifs**2)
        elif classname.find('Linear') != -1:
            m.weight.data.normal_(0.0, self.sigma_net**2)
            m.bias.data.fill_(0.00001)

    def forward(self, x):
        x = self.conv_1(x)
        x = self.nonlin(x)
        x = self.pool_1(x)
        # Remove channel dimension
        x = x.view(-1, self.num_motifs)
        x = self.fc(x)
        return x

    def training_step(self, batch, batch_idx):
        seqs, vals = batch
        outputs = self.forward(seqs)
        loss = self.loss(outputs, vals)
        self.log("ptl/train_loss", loss)
        return {'loss': loss}

    def validation_step(self, batch, batch_idx):
        seqs, vals = batch
        outputs = self.forward(seqs)
        loss = self.loss(outputs, vals)
        return {"val_loss": loss}

    def training_epoch_end(self, outputs):
        avg_loss = torch.stack([x["loss"] for x in outputs]).mean()
        self.log("ptl/epoch_train_loss", avg_loss)

    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack([x["val_loss"] for x in outputs]).mean()
        self.log("ptl/val_loss", avg_loss)

    def test_step(self, batch, batch_idx):
        seqs, vals = batch
        outputs = self.forward(seqs)
        loss = self.loss(outputs, vals)
        self.log('test_loss', loss)

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        return optimizer
