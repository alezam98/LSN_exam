#!/usr/bin/env python
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

seed = 0
np.random.seed(seed)
tf.random.set_seed(seed)

NDATA = 10000
SIGMA = 1


def f(x):
    return 2.*x + 1
    

def generate_data():
    x_train = np.random.uniform(-1, 1, NDATA)
    x_test = np.random.uniform(-1, 1, int(NDATA/5))

    y_train = np.random.normal(f(x_train), SIGMA, size=NDATA)
    y_test = np.random.normal(f(x_test), SIGMA, size=int(NDATA/5))
    y_target = f(x_test)

    return x_train, y_train, x_test, y_test, y_target


def generate_model():
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(units=1, input_shape=(1,), activation="linear"))

    model.compile(
        loss = tf.keras.losses.MeanSquaredError(), 
        optimizer = tf.keras.optimizers.SGD()
    )

    return model


def plot_results(model, x_test, y_test, history, metric, filename="results.png", show_bool=False, save_bool=False):
    x = np.linspace(-1.5, 1.5, 100)
    predictions = model.predict(x)

    train_metrics = history.history[metric]
    val_metrics = history.history['val_'+metric]
    epochs = range(0, len(train_metrics))
    epoch = np.where(val_metrics == np.min(val_metrics))[0][0]

    plt.rcParams.update({'font.size':13})
    fig, (ax1, ax2) = plt.subplots(figsize=(16, 5), ncols=2)
    
    ax = ax1
    ax.set_title("Model predictions")
    ax.scatter(x_test, y_test, color="green", marker=".", alpha=0.5, label="test data")
    ax.scatter(x, predictions, color="blue", label="predictions")
    ax.plot(x, f(x), color="red", label="target")
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.legend()

    ax = ax2
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



def print_weights(model, filename="weights.dat", save_bool=False):
    print("\nModel parameters:")
    print(f'm = {model.get_weights()[0][0][0]}, c = {model.get_weights()[1][0]}')
    print()

    if save_bool:
        with open(filename, 'w') as fp:
            fp.write(f"{model.get_weights()[0][0][0]}\n")
            fp.write(f"{model.get_weights()[1][0]}")



def main():

    # data loading
    x_train, y_train, x_test, y_test, y_target = generate_data()

    print("\nTrain shape:")
    print(f"x: {x_train.shape}, y: {y_train.shape}")
    print("Test shape:")
    print(f"x: {x_test.shape}, y: {y_test.shape}\n")

    # model generation
    model = generate_model()
    
    # fit
    history = model.fit(
        x_train, y_train,
        batch_size = 32,
        validation_split = 0.4,
        epochs = 50,
        shuffle = True,
        callbacks = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)]
    )
    metric = model.evaluate(x_test, y_target, verbose=1)

    # results
    plot_results(model, x_test, y_test, history, "loss", "../data/results_01.png", show_bool=False, save_bool=True)
    print_weights(model, filename="../data/weights_01.dat", save_bool=True)
    print_evaluation(metric, filename="../data/test_metric_01.dat", save_bool=True)
          


if __name__ == "__main__":
    main()