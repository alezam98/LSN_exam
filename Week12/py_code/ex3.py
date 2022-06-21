#!/usr/bin/env python
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from numpy import asarray
from os import listdir
from os.path import isfile, isdir, join
import re

seed = 0
np.random.seed(seed)
tf.random.set_seed(seed)

DNN_BEST = 'adam'
LOAD = False


def load_models():
    DNN_model = tf.keras.models.load_model('../models/DNN_'+DNN_BEST+'_model')
    CNN_model = tf.keras.models.load_model('../models/CNN_model')

    models = {
        "DNN": DNN_model,
        "CNN": CNN_model
    }
    return models


def load_images():
    for type_ in ['normal', 'weird']:
        mypath = '../new_images/'+type_
        files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
        files = [re.split('\.',f)[0] for f in files]

        images = np.zeros(shape=(len(files), 28, 28))
        labels = np.zeros(shape=(len(files),))

        for i in range(len(files)):
            image = asarray(Image.open('../new_images/'+type_+'/'+files[i]+'.png'))
            images[i] = image
            labels[i] = int(files[i])

        if type_ == 'normal': normal_images, normal_labels = images, labels
        elif type_ == 'weird': weird_images, weird_labels = images, labels

    return normal_images, normal_labels, weird_images, weird_labels


def generate_DA_model():
    model = tf.keras.models.Sequential()

    # Data Augmentation layers
    model.add(tf.keras.layers.RandomRotation(0.3, input_shape=(28, 28, 1)))
    model.add(tf.keras.layers.RandomZoom(0.2))
    model.add(tf.keras.layers.RandomTranslation(height_factor=(-0.4, 0.4), width_factor=(-0.4, 0.4)))

    model.add(tf.keras.layers.Rescaling(1./255, input_shape=(28, 28, 1)))
    model.add(tf.keras.layers.Conv2D(32, (3, 3), activation='relu'))
    model.add(tf.keras.layers.MaxPool2D((2, 2)))
    model.add(tf.keras.layers.Conv2D(16, (3, 3), activation='relu'))
    model.add(tf.keras.layers.MaxPool2D((2, 2)))

    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(units=64, activation="relu"))
    model.add(tf.keras.layers.Dense(units=32, activation="relu"))
    model.add(tf.keras.layers.Dense(units=10, activation="softmax"))

    model.compile(
        optimizer = 'adam',
        loss = 'sparse_categorical_crossentropy',
        metrics = ['accuracy']
    )

    return model


def train(model, train_images, train_labels):
    history = model.fit(
        train_images, train_labels,
        validation_split = 0.3,
        batch_size = 32,
        epochs = 100,
        callbacks = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=20)]
    )
    return history


def plot_predictions(models, test_images, test_labels, filename='predictions.png', save_bool=False, show_bool=False):
    predictions = {}
    for model in models: predictions[model] = models[model].predict(test_images)

    plt.rcParams.update({'font.size':7})
    fig = plt.figure(figsize=(8,5))
    for i in range(10):
        ax = plt.subplot(2, 5, i+1)
        
        ax.imshow(test_images[i], cmap='gray')
        if len(models) == 2:
            ax.set_title(f"Number: {int(test_labels[i])}\nDNN prediction: {np.argmax(predictions['DNN'][i])}\nCNN prediction: {np.argmax(predictions['CNN'][i])}")
        elif len(models) == 3:
            ax.set_title(f"Number: {int(test_labels[i])}\nDNN prediction: {np.argmax(predictions['DNN'][i])}\nCNN prediction: {np.argmax(predictions['CNN'][i])}\nDA prediction: {np.argmax(predictions['DA'][i])}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    
    fig.tight_layout()
    if save_bool: plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    if show_bool: plt.show()


def plot_metrics(history, metrics, filename='metrics.png', save_bool=False, show_bool=False):
    num_plots = len(metrics)

    plt.rcParams.update({'font.size':10})
    fig = plt.figure(figsize=(8*num_plots, 5))
    for metric in metrics:
        train_metric = history.history[metric]
        valid_metric = history.history['val_'+metric]
        epochs = np.arange(len(valid_metric))

        ax = plt.subplot(1, num_plots, metrics.index(metric)+1)
        ax.plot(epochs, train_metric, color='red', label='Training')
        ax.plot(epochs, valid_metric, color='blue', label='Validation')

        if metric == 'loss': epoch = np.where(valid_metric == np.min(valid_metric))[0][0]
        ax.scatter(epoch, valid_metric[epoch], color='blue', label=f'Best model (epoch={epoch})')

        ax.set_title('Training and validation '+metric)
        ax.set_xlabel('epochs')
        ax.set_ylabel(metric)
        ax.legend()

    fig.tight_layout()
    if save_bool: plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    if show_bool: plt.show()


def print_evaluation(metrics, dir_=None, save_bool=False):

    print("\DA model evaluated on the test set:")
    print(f"loss={metrics[0]}, accuracy={metrics[1]}")
    print()

    if save_bool:
        with open(dir_+'test_metrics.dat', 'w') as fp:
            fp.write(f'{metrics[0]}\n')
            fp.write(f'{metrics[1]}\n')



def main():

    # load old models and new images
    models = load_models()
    normal_images, normal_labels, weird_images, weird_labels = load_images()
    dir_ = '../data/ex3/'

    # test of the old models on the normal images
    plot_predictions(models, normal_images, normal_labels, filename=dir_+'first_normal_prediction.png', save_bool=True, show_bool=True)
    plot_predictions(models, weird_images, weird_labels, filename=dir_+'first_weird_prediction.png', save_bool=True, show_bool=True)


    # generating new CNN model with data augmentation (DA_model)
    mypath = '../models'
    files = [f for f in listdir(mypath) if isdir(join(mypath, f))]
    mask = np.any([f == 'DA_model' for f in files])
    
    if np.all([mask, LOAD]):
        DA_model = tf.keras.models.load_model('../models/DA_model') 
        models["DA"] = DA_model
    else:
        (train_images, train_labels), (test_images, test_labels) = tf.keras.datasets.mnist.load_data()

        print("\nTrain shape:")
        print(f"images: {train_images.shape}, labels: {train_labels.shape}")
        print("Test shape:")
        print(f"images: {test_images.shape}, labels: {test_labels.shape}\n")
        
        DA_model = generate_DA_model()
        history = train(DA_model, train_images, train_labels)
        metrics = DA_model.evaluate(test_images, test_labels, verbose=2)
        DA_model.save('../models/DA_model')
        models["DA"] = DA_model
    
        plot_metrics(history, ['loss', 'accuracy'], filename=dir_+'metrics.png', save_bool=True, show_bool=True)
        print_evaluation(metrics, dir_=dir_, save_bool=True)


    # test of all the models on both normal and weird images
    plot_predictions(models, normal_images, normal_labels, filename=dir_+'second_normal_prediction.png', save_bool=True, show_bool=True)
    plot_predictions(models, weird_images, weird_labels, filename=dir_+'second_weird_prediction.png', save_bool=True, show_bool=True)




if __name__ == '__main__':
    main()
