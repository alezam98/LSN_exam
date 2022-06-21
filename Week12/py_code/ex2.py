#!/usr/bin/env python
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

seed = 0
np.random.seed(seed)
tf.random.set_seed(seed)


def generate_CNN_model():
    model = tf.keras.models.Sequential()

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



def plot_predictions(model, test_images, test_labels, filename='predictions.png', save_bool=False, show_bool=False):
    np.random.seed(seed+1)
    predictions = model.predict(test_images)

    plt.rcParams.update({'font.size':7})
    fig = plt.figure(figsize=(5,5))
    for i in range(9):
        ax = plt.subplot(3, 3, i+1)
        
        iim = np.random.randint(0, len(test_images))
        ax.imshow(test_images[iim], cmap='gray')
        ax.set_title(f"Number: {test_labels[iim]}\nPrediction: {np.argmax(predictions[iim])}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    
    fig.tight_layout()
    if save_bool: plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    if show_bool: plt.show()
    np.random.seed(seed)




def plot_metrics(history, metrics, filename='metrics.png', save_bool=False, show_bool=False):
    num_plots = len(metrics)

    plt.rcParams.update({'font.size':13})
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

    print("\nCNN model evaluated on the test set:")
    print(f"loss={metrics[0]}, accuracy={metrics[1]}")
    print()

    if save_bool:
        with open(dir_+'test_metrics.dat', 'w') as fp:
            fp.write(f'{metrics[0]}\n')
            fp.write(f'{metrics[1]}\n')





def main():

    # loading data
    (train_images, train_labels), (test_images, test_labels) = tf.keras.datasets.mnist.load_data()
    
    print("\nTrain shape:")
    print(f"images: {train_images.shape}, labels: {train_labels.shape}")
    print("Test shape:")
    print(f"images: {test_images.shape}, labels: {test_labels.shape}\n")

    ## exercise 12.02
    CNN_model = generate_CNN_model()
    history = train(CNN_model, train_images, train_labels)
    metrics = CNN_model.evaluate(test_images, test_labels, verbose=2)
    CNN_model.save('../models/CNN_model')

    # results
    dir_ = '../data/ex2/' 
    plot_metrics(history, ['loss', 'accuracy'], filename=dir_+'metrics.png', save_bool=True, show_bool=False)
    plot_predictions(CNN_model, test_images, test_labels, filename=dir_+'prediction.png', save_bool=True, show_bool=False)
    print_evaluation(metrics, dir_=dir_, save_bool=True)
    


if __name__ == "__main__":
    main()