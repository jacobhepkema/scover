import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np

## FUNCTIONS ================================================
def loguniform(a=0, b=1, size=None):
    X = 10 ** ((np.log10(b) - np.log10(a)) * 
               np.random.random_sample(size) + 
               np.log10(a))
    return X

def initialize_learning_rate(epsilon_min, epsilon_max):
    learning_rate = loguniform(a = epsilon_min, b = epsilon_max)
    return learning_rate

def get_errtest_len(num_learningsteps, num_errtest):
    errtest_len = 0
    for i in range(num_learningsteps):
        if (i+1) % num_errtest == 0:
            errtest_len += 1
    return errtest_len

## NEURAL NETWORK CLASS DEFINITION ==========================
class ConvNet(nn.Module):
    # Initialization method
    def __init__(self,  
                 optimizer_name,
                 criterion_name, 
                 epsilon_min, 
                 epsilon_max, 
                 sigma_motifs_min, 
                 sigma_motifs_max, 
                 sigma_net_min, 
                 sigma_net_max, 
                 d, 
                 m, 
                 n, 
                 cells):    
        super(ConvNet, self).__init__()
        
        # Weight initialization method
        def weights_init(m):
            classname = m.__class__.__name__
            if classname.find('Conv2d') != -1:
                self.sigma_motifs = loguniform(sigma_motifs_min, sigma_motifs_max)
                m.weight.data.normal_(0.0, self.sigma_motifs**2)
            elif classname.find('Linear') != -1:
                self.sigma_net = loguniform(sigma_net_min, sigma_net_max)
                m.weight.data.normal_(0.0, self.sigma_net**2)
                m.bias.data.fill_(0.00001)
        
        # First convolution and pooling
        self.conv_1 = nn.Conv2d(1, d, (m, 4), bias = False)
        self.pool_1 = nn.MaxPool2d(((n - m + 1), 1))
        
        # Fully-connected layer
        self.fc = nn.Linear(d, cells, bias = True)
        
        self.d = d
        self.m = m
        self.n = n
        self.cells = cells
        
        self.beta = [loguniform(1e-15,1e-3),  
                     loguniform(1e-10,1e-3)]
        
        # Initialize network weights
        self.learning_rate = initialize_learning_rate(epsilon_min, epsilon_max)
        optimizers = {
            "Adagrad": optim.Adagrad([{'params': self.conv_1.parameters(),
                                      'weight_decay': self.beta[0]},
                                      {'params': self.fc.parameters(),
                                      'weight_decay': self.beta[1]}],
                                      lr=self.learning_rate),
            "Adam": optim.Adam([{'params': self.conv_1.parameters(),
                'weight_decay': self.beta[0]},
                {'params': self.fc.parameters(),
                    'weight_decay': self.beta[1]}],
                lr=self.learning_rate),
            "SGD": optim.SGD([{'params': self.conv_1.parameters(), 
                              'weight_decay': self.beta[0]}, 
                              {'params': self.fc.parameters(),
                              'weight_decay': self.beta[1]}],
                             lr = self.learning_rate,
                             momentum=0.97)
        }
        self.optimizer = optimizers.get(optimizer_name, "Invalid 'optimizer_name'")
        criteria = {
            "MSE": nn.MSELoss().double(),
            "L1": nn.SmoothL1Loss().double()
        }
        self.criterion = criteria.get(criterion_name, "Invalid 'criterion_name'")
        self.apply(weights_init)
        
    def forward(self, x):
        x = self.conv_1(x)
        x = F.relu(x)
        x = self.pool_1(x)
        # Remove dimension of channels
        x = x.view(-1, self.d)
        x = self.fc(x)
        return x
    
## NEURAL NETWORK CLASS DEFINITION ==========================
class BestInitialConvNet(nn.Module):
    # Initialization method
    def __init__(self,  
                 optimizer_name,
                 criterion_name, 
                 learning_rate, 
                 sigma_motifs, 
                 sigma_net,
                 d, 
                 m, 
                 n, 
                 cells):    
        super(BestInitialConvNet, self).__init__()
        
        self.learning_rate = learning_rate
        self.sigma_motifs = sigma_motifs
        self.sigma_net = sigma_net

        # Weight initialization method
        def weights_init(m):
            classname = m.__class__.__name__
            if classname.find('Conv2d') != -1:
                m.weight.data.normal_(0.0, self.sigma_motifs**2)
            elif classname.find('Linear') != -1:
                m.weight.data.normal_(0.0, self.sigma_net**2)
                m.bias.data.fill_(0.00001)
        
        # First convolution and pooling
        self.conv_1 = nn.Conv2d(1, d, (m, 4), bias = False)
        self.pool_1 = nn.MaxPool2d(((n - m + 1), 1))
        
        # Fully-connected layer
        self.fc = nn.Linear(d, cells, bias = True)
        
        self.d = d
        self.m = m
        self.n = n
        self.cells = cells
        
        self.beta = [loguniform(1e-15,1e-3),  
                     loguniform(1e-10,1e-3)]
        
        # Initialize network weights
        optimizers = {
            "Adagrad": optim.Adagrad([{'params': self.conv_1.parameters(),
                                      'weight_decay': self.beta[0]},
                                      {'params': self.fc.parameters(),
                                      'weight_decay': self.beta[1]}],
                                      lr=self.learning_rate),
            "Adam": optim.Adam([{'params': self.conv_1.parameters(),
                'weight_decay': self.beta[0]},
                {'params': self.fc.parameters(),
                    'weight_decay': self.beta[1]}],
                lr=self.learning_rate),
            "SGD": optim.SGD([{'params': self.conv_1.parameters(), 
                              'weight_decay': self.beta[0]}, 
                              {'params': self.fc.parameters(),
                              'weight_decay': self.beta[1]}],
                             lr = self.learning_rate,
                             momentum=0.97)
        }
        self.optimizer = optimizers.get(optimizer_name, "Invalid 'optimizer_name'")
        criteria = {
            "MSE": nn.MSELoss().double(),
            "L1": nn.SmoothL1Loss().double()
        }
        self.criterion = criteria.get(criterion_name, "Invalid 'criterion_name'")
        self.apply(weights_init)
        
    def forward(self, x):
        x = self.conv_1(x)
        x = F.relu(x)
        x = self.pool_1(x)
        # Remove dimension of channels
        x = x.view(-1, self.d)
        x = self.fc(x)
        return x
