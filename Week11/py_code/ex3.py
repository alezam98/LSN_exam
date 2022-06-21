#!/usr/bin/env python
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

seed = 0
np.random.seed(seed)
tf.random.set_seed(seed)

NDATA = 100000
SIGMA = 1


def f(x1, x2):
    return np.sin(x1**2. + x2**2.)


def generate_data():
    x_train = np.random.uniform(-1.5, 1.5, (NDATA, 2))
    x_test = np.random.uniform(-1.5, 1.5, (int(NDATA/5), 2))

    y_train = np.random.normal(f(x_train[:, 0], x_train[:, 1]), SIGMA, size=NDATA)
    y_target = f(x_test[:, 0], x_test[:, 1])
    
    return x_train, y_train, x_test, y_target


def generate_model():
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(units=64, activation="relu", input_dim=2))
    model.add(tf.keras.layers.Dense(units=32, activation="relu"))
    model.add(tf.keras.layers.Dense(units=1, activation="linear"))

    model.compile(
        loss = tf.keras.losses.MeanSquaredError(), 
        optimizer = tf.keras.optimizers.SGD(0.1)
    )

    return model


def plot_results(model, history, metric, filename="results.png", show_bool=False, save_bool=False):
    x = np.random.uniform(-1.5, 1.5, (100, 2))
    predictions = model.predict(x)

    X, Y = np.meshgrid(np.linspace(-1.5, 1.5, 100), np.linspace(-1.5, 1.5, 100))
    Z = f(X, Y)

    train_metrics = history.history[metric]
    val_metrics = history.history['val_'+metric]
    epochs = range(0, len(train_metrics))
    epoch = np.where(val_metrics == np.min(val_metrics))[0][0]

    plt.rcParams.update({'font.size':13})
    fig = plt.figure(figsize=(16, 5))
    
    ax = plt.subplot(1, 2, 1, projection="3d")
    ax.set_title("Model predictions")
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.scatter(x[:, 0], x[:, 1], predictions, marker="o", color="red")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("f(x,y)")
    #ax.view_init(60, 35)

    ax = plt.subplot(1, 2, 2)
    ax.plot(epochs, train_metrics, linestyle="-", color="red", label="Training")
    ax.plot(epochs, val_metrics, linestyle="-", color="blue", label="Validation")
    ax.scatter(epoch, val_metrics[epoch], color='blue', label=f'Best model (epoch={epoch})')
    ax.set_title('Training and validation '+ metric)
    ax.set_xlabel("epochs")
    ax.set_ylabel(metric)
    ax.legend()

    fig.tight_layout()
    if save_bool: plt.savefig(filename, bbox_inches='tight', pad_inches=0.05)
    if show_bool: plt.show()



def print_evaluation(metric, filename="test_metric.dat", save_bool=False):

    print("\nModel evaluated on the test set:")
    print(f"loss={metric}")
    print()

    if save_bool:
        with open(filename, 'w') as fp:
            fp.write(f'{metric}')



def main():

    # data loading
    x_train, y_train, x_test, y_target = generate_data()

    print("\nTrain shape:")
    print(f"x: {x_train.shape}, y: {y_train.shape}\n")

    # model generation
    model = generate_model()
    
    # fit
    history = model.fit(
        x_train, y_train,
        batch_size = 32,
        validation_split = 0.4,
        epochs = 50,
        shuffle = True,
        callbacks = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=20)]
    )
    metric = model.evaluate(x_test, y_target, verbose=1)

    # results
    plot_results(model, history, "loss", "../data/results_03.png", show_bool=False, save_bool=True)
    print_evaluation(metric, filename="../data/test_metric_03.dat", save_bool=True)


if __name__ == "__main__":
    main()