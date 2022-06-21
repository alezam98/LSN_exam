#!/usr/bin/env python
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

seed = 0
np.random.seed(seed)
tf.random.set_seed(seed)


def generate_DNN_model(optimizer):
    model = tf.keras.models.Sequential()

    model.add(tf.keras.layers.Rescaling(1./255., input_shape=(28, 28)))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(units=64, activation="relu"))
    model.add(tf.keras.layers.Dense(units=32, activation="relu"))
    model.add(tf.keras.layers.Dense(units=10, activation="softmax"))

    model.compile(
        optimizer = optimizer,
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



def compare_metrics(histories, metrics, dir_=None, save_bool=False, show_bool=False):
    best_epoch = {}

    plt.rcParams.update({'font.size':13})
    for metric in metrics:
        fig = plt.figure(figsize=(16, 5))

        if metric == 'loss': fig.suptitle("Loss function, DNN models", y=0.95)
        elif metric == 'accuracy': fig.suptitle("Accuracy, DNN models", y=0.95)

        ax = plt.subplot(1, 2, 1)
        for optimizer in histories:
            train_metric = histories[optimizer].history[metric]
            epochs = np.arange(len(train_metric))
            ax.plot(epochs, train_metric, label=optimizer)

        ax.set_title('Training set')
        ax.set_xlabel('epochs')
        ax.set_ylabel(metric)
        ax.legend()


        ax = plt.subplot(1, 2, 2)
        for optimizer in histories:
            valid_metric = histories[optimizer].history['val_'+metric]
            epochs = np.arange(len(valid_metric))
            ax.plot(epochs, valid_metric, label=optimizer)
            
            if metric == 'loss': best_epoch[optimizer] = np.where(valid_metric == np.min(valid_metric))[0][0]
            ax.scatter(best_epoch[optimizer], valid_metric[best_epoch[optimizer]], label=f'{optimizer}, best (epoch={best_epoch[optimizer]})')

        ax.set_title('Validation set')
        ax.set_xlabel('epochs')
        ax.set_ylabel(metric)
        ax.legend()

        
        fig.tight_layout()
        if save_bool: plt.savefig(dir_+metric+'_metrics.png', bbox_inches='tight', pad_inches=0.05)
        if show_bool: plt.show()



def print_evaluation(metrics, dir_=None, save_bool=False):

    print("\nDNN models evaluated on the test set:")
    for optimizer in metrics:
        print(f"{optimizer}: loss={metrics[optimizer][0]}, accuracy={metrics[optimizer][1]}")
    print()

    if save_bool:
        with open(dir_+'test_metrics.dat', 'w') as fp:
            for optimizer in metrics:
                fp.write(f'{optimizer}\n')
                fp.write(f'{metrics[optimizer][0]}\n')
                fp.write(f'{metrics[optimizer][1]}\n')





def main():

    # loading data
    (train_images, train_labels), (test_images, test_labels) = tf.keras.datasets.mnist.load_data()
    
    print("\nTrain shape:")
    print(f"images: {train_images.shape}, labels: {train_labels.shape}")
    print("Test shape:")
    print(f"images: {test_images.shape}, labels: {test_labels.shape}\n")

    ## exercise 12.01
    optimizers = ['adam', 'sgd', 'adamax']
    metrics = {}
    histories = {}
    dir_ = '../data/ex1/'

    for optimizer in optimizers:
        print(f'\n\nOPTIMIZER: {optimizer}\n')
        DNN_model = generate_DNN_model(optimizer)
        history = train(DNN_model, train_images, train_labels)
        plot_metrics(history, ['loss', 'accuracy'], filename=dir_+optimizer+'_metrics.png', save_bool=True, show_bool=False)
        plot_predictions(DNN_model, test_images, test_labels, filename=dir_+optimizer+'_prediction.png', save_bool=True, show_bool=False)
        
        metrics[optimizer] = DNN_model.evaluate(test_images, test_labels, verbose=2)
        histories[optimizer] = history
        DNN_model.save(f'../models/DNN_{optimizer}_model')

    # results
    compare_metrics(histories, ['loss', 'accuracy'], dir_=dir_, save_bool=True, show_bool=False)
    print_evaluation(metrics, dir_=dir_, save_bool=True)
    


if __name__ == "__main__":
    main()