from keras.layers import Dense, Dropout  # Conv2D, MaxPooling2D, UpSampling2D,
from keras import Input, Model
import numpy as np
import pandas as pd
import time
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import GridSearchCV
from matplotlib import pyplot as plt

example_data = pd.read_table(
    filepath_or_buffer="S:/KJP_Biolabor/Projects/RAISE_GENIC/data/train_SNPs_pathway_1.raw", sep=" ")
x_train = example_data.drop(["FID", "IID", "SEX", "MAT", "PAT", "PHENOTYPE"], axis=1)


# check handling of missingness [imputed data wont have missingness
x_train.fillna(0, inplace=True)

# normalize between 0 and 1
x_train = x_train / 2

snplength = len(x_train.columns)
minlayer = 0.05
patients = len(x_train)

activation = ['relu', "sigmoid"]
activation_last = ['sigmoid']
dropout_rate = [0, 0.2, 0.5]
layers = [0, 1, 2] # = layers between input & bneck, symmetrically
optimizer = ['SGD', 'Adam'] # check adam properties
loss = ['mean_squared_error']
batches = [20]
epochs = [50]

def create_model(activation,
                 dropout_rate,
                 layers,
                 activation_last,
                 optimizer,
                 loss,
                 vectorlen,
                 **args):
    input_snps = Input(shape=(vectorlen,))
    encoded = Dropout(dropout_rate)(input_snps)

    if layers > 0:
        encoded = Dense(vectorlen / 2, activation=activation)(encoded)
    if layers == 2:
        encoded = Dense(vectorlen / 4, activation=activation)(encoded)

    decoded = Dense(int(np.ceil(vectorlen * minlayer)), activation=activation)(encoded)

    if layers == 2:
        decoded = Dense(vectorlen / 4, activation=activation)(decoded)
    if layers != 0:
        decoded = Dense(vectorlen / 2, activation=activation)(decoded)

    decoded = Dense(vectorlen, activation=activation_last)(decoded)
    autoencoder = Model(input_snps, decoded)
    autoencoder.compile(optimizer=optimizer, loss=loss)
    return autoencoder


if __name__ == "__main__":
    param_grid = dict(optimizer=optimizer,
                      dropout_rate=dropout_rate,
                      activation=activation,
                      activation_last=activation_last,
                      loss=loss,
                      layers=layers,
                      batch_size=batches,
                      epochs=epochs,
                      vectorlen=[snplength])

    # model = create_model(optimizer=optimizer[0], dropout_rate=dropout_rate[0],
    #                      activation=activation[0], activation_last=activation_last[0],
    #                      loss=loss[0],
    #                      layers=layers[0],
    #                      vectorlen=snplength)
    # model.fit(x_train, x_train, epochs=epochs[0], batch_size=batches[0])


    ## new KerasRegressor version has problems handling variables stick to deprecated for now
    model = KerasRegressor(build_fn=create_model)

    grid = GridSearchCV(estimator=model, param_grid=param_grid, n_jobs=10, cv=3, verbose=1)

    start_time = time.time()
    grid.fit(x_train, x_train)

    print("--- %s seconds ---" % (time.time() - start_time))
    print("--- %s models ---" % len(grid.cv_results_["mean_test_score"]))
    print("--- %s SNPs ---" % snplength)
    print("--- %s patients ---" % patients)

    #plot learing over epochs

    grid.best_score_
    bestmodel = grid.best_estimator_

    model = create_model(**grid.best_params_)
    hist = model.fit(x_train, x_train, batch_size=grid.best_params_["batch_size"],
                     epochs=grid.best_params_["epochs"])

    plt.plot(hist.history["loss"])


if False:
    class Hypertrainer:
        def __init__(self, param_grid, x_train):
            param_grid=param_grid
            x_train=x_train

        def create_model(self,
                         activation,
                         dropout_rate,
                         layers,
                         activation_last,
                         optimizer,
                         loss,
                         vectorlen,
                         **args):
            input_snps = Input(shape=(vectorlen,))
            encoded = Dropout(dropout_rate)(input_snps)

            if layers > 0:
                encoded = Dense(snplength / 2, activation=activation)(encoded)
            if layers == 2:
                encoded = Dense(snplength / 4, activation=activation)(encoded)
            decoded = Dense(int(np.ceil(snplength * minlayer)), activation=activation)(encoded)

            if layers == 2:
                decoded = Dense(snplength / 4, activation=activation)(decoded)
            if layers != 0:
                decoded = Dense(snplength / 2, activation=activation)(decoded)

            decoded = Dense(snplength, activation=activation_last)(decoded)
            autoencoder = Model(input_snps, decoded)
            autoencoder.compile(optimizer=optimizer, loss=loss)
            return autoencoder



