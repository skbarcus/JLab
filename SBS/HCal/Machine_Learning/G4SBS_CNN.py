import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

from matplotlib.ticker import MaxNLocator #Able to force plots to use integer ticks only. 

#tf.logging.set_verbosity(tf.logging.ERROR)#Suppress Tensorflow warnings.

def plot_learning_curve(history):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(history['loss'], label='training')
    ax[0].plot(history['val_loss'], label='validation')
    ax[0].set_title("Model losses")
    ax[0].set_xlabel("Epoch")
    ax[0].set_ylabel("Loss")
    ax[0].set_yscale('log')
    ax[0].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[0].legend()
    
    ax[1].plot(history['accuracy'], label='training')
    ax[1].plot(history['val_accuracy'], label='validation')
    ax[1].set_title("Model accuracy")
    ax[1].set_xlabel("Epoch")
    ax[1].set_ylabel("Accuracy")
    ax[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[1].legend()
    plt.show()

data = "/home/skbarcus/JLab/SBS/HCal/Analysis/Simulation/G4SBS/rootfiles/"

fnucl = np.load(data+"g4sbs_kin7_fnucl.npy")
front_nhits = np.load(data+"g4sbs_kin7_front_nhits.npy")
nhits = np.load(data+"g4sbs_kin7_nhits.npy")
pmt = np.load(data+"g4sbs_kin7_pmt.npy", allow_pickle=True)
pmt_edep = np.load(data+"g4sbs_kin7_pmt_edep.npy", allow_pickle=True)
tavg = np.load(data+"g4sbs_kin7_tavg.npy", allow_pickle=True)

