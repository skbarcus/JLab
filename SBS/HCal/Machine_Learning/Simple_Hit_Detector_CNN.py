import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

#tf.logging.set_verbosity(tf.logging.ERROR)#Suppress Tensorflow warnings.

def plot_learning_curve(history):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(history['loss'], label='training')
    ax[0].plot(history['val_loss'], label='validation')
    ax[0].set_title("Model losses")
    ax[0].set_xlabel("Epoch")
    ax[0].set_ylabel("Loss")
    ax[0].set_yscale('log')
    ax[0].legend()
    
    ax[1].plot(history['accuracy'], label='training')
    ax[1].plot(history['val_accuracy'], label='validation')
    ax[1].set_title("Model accuracy")
    ax[1].set_xlabel("Epoch")
    ax[1].set_ylabel("Accuracy")
    ax[1].legend()
    plt.show()
    # plt.plot(history["loss"], label="training loss")
    # plt.plot(history["val_loss"], label="validation loss")
    # plt.xlabel("Epoch")
    # plt.ylabel("Loss")
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

adc = np.load("/home/skbarcus/Machine_Learning/Tutorials/Data/adc_integrals_run820.npy")
hits = np.load("/home/skbarcus/Machine_Learning/Tutorials/Data/hits_run820.npy")

print('adc matrix = ', adc)
print('adc.shape = ', adc.shape)
adc = np.reshape(adc, (adc.shape[0], 12, 12, 1))
print('adc.shape = ', adc.shape)
print('hits matrix = ', hits)
print('hits.shape = ', hits.shape)

#Split data into training and test.
adc_train, adc_test, hits_train, hits_test = train_test_split(adc, hits, test_size=0.25)

print('adc_train.shape, adc_test.shape = ', adc_train.shape, adc_test.shape)
print('hits_train.shape, hits_test.shape = ', hits_train.shape, hits_test.shape)

unique, counts = np.unique(hits_train, return_counts=True)
counts = dict(zip(unique, counts))
print('Counts = ', counts)

model = tf.keras.Sequential() #Define the model object

model.add(tf.keras.layers.Conv2D(filters=32,kernel_size=3,activation='relu',input_shape=adc.shape[1:]))#Shape of a single image.
model.add(tf.keras.layers.MaxPool2D(2,2))
model.add(tf.keras.layers.Flatten())
model.add(tf.keras.layers.Dense(128, activation='relu'))
model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
model.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy']
)
"""
#Dense NN (not CNN layers).
model.add(tf.keras.layers.Dense(128, input_shape=(144,), activation="relu")) #Add the hidden layer
model.add(tf.keras.layers.Dense(64, activation="relu"))
#model.add(tf.keras.layers.Dense(32, activation="relu"))
model.add(tf.keras.layers.Dense(1))
"""

#model.compile(tf.keras.optimizers.Adam(lr=0.01),loss=tf.keras.losses.CategoricalCrossentropy()) #Adam optimizer and mean squared error loss
#model.compile(tf.keras.optimizers.Adam(lr=0.0001),loss=tf.keras.losses.MeanSquaredError(),metrics=['accuracy'])

#Print summary of model.
model.summary()

results = model.fit(adc_train, hits_train, epochs=20, batch_size=64, validation_split=0.3333)

#prediction = model.predict()
scores = model.evaluate(adc_test, hits_test)

print("test loss, test acc:", scores)

plot_learning_curve(results.history)
