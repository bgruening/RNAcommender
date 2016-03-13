#!/usr/bin/env python

from __future__ import print_function
import argparse
import cPickle

from data import TrainDataset
from model import Model

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

class Predictor():
    def __init__(self,predict_dataset,trained_model,output=None,verbose=True):
        """
        Params
        ------
        predict_dataset : data.PredictDataset
            Dataset containing the examples to predict.

        trained_model : str
            File name of the trained model.

        output : str (default : None)
            Output file. If None then STDOUT

        verbose : bool (default : True)
            Print information at STDOUT.
        """
        self.predict_dataset = predict_dataset
        f = open(trained_model)
        self.model = cPickle.load(f)
        f.close()
        self.output = output
        self.verbose = verbose

    def predict():
        "Predict interaction values"
        if self.verbose:
            print("Predicting interactions...", end=' ')
            sys.stdout.flush()
        # predict the y_hat
        (P,P_names,R,R_names) = self.predict_dataset
        y_hat = model.predict(P,R)
        # sort the interactions according to y_hat
        ordering = sorted(range(len(y_hat)),key=lambda x:y_hat[x],reverse=True)
        P_names = P_names[ordering]
        R_names = R_names[ordering]
        y_hat = y_hat[ordering]

        if self.verbose:
            print("Done.")
            sys.stdout.flush()

        # output to STDOUT
        if self.output is None:
            print("RBP\ttarget\ty_hat")
            for (p,r,s) in izip(P_names,R_names,y_hat):
                print("%s\t%s\t%.3f" % (p,r,s))
        # output to file
        else:
            nf = open(self.output,"w")
            nf.write("RBP\ttarget\ty_hat\n")
            for (p,r,s) in izip(P_names,R_names,y_hat):
                nf.write("%s\t%s\t%.3f\n" % (p,r,s))
            nf.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('Fp', metavar='Fp', type=str,
                        help="""HDF5 file with the protein features.""")
    parser.add_argument('Fr', metavar='Fr', type=str,
                        help="""HDF5 file with the RNA features.""")
    parser.add_argument('model', metavar='model', type=str,
                        help="""Trained model to use for the prediction.""")
    parser.add_argument('--to-predict', metavar='to_predict', nargs='+', required=False, default=None,
                        help="""Space separated list of proteins from Fp to predict.""")
    parser.add_argument('--standardize-Fp', dest='standardize_Fp', action='store_true', default=False,
                        help="""Standardize protein features.""")
    parser.add_argument('--standardize-Fr', dest='standardize_Fr', action='store_true', default=False,
                        help="""Standardize RNA features.""")
    parser.add_argument('--output', metavar='output', type=str, default=None,
                        help="""File name to save predictions. Default is STDOUT.""")
    parser.add_argument('--quiet', dest='quiet', action='store_true', default=False,
                        help="""Silence STDOUT.""")

    args = parser.parse_args()

    def feature_size(store_name):
        """Number of features"""
        store = pd.io.pytables.HDFStore(store_name)
        a = store.features
        store.close()
        return a.shape[0]

    # Define and instantiate dataset
    D = PredictDataset(Fp=args.Fp,Fr=args.Fr,to_predict=args.to_predict,
        standardize_proteins=args.standardize_Fp, standardize_rnas=args.standardize_Fr,
        verbose=(not args.quiet))
    dataset = D.load()

    # Define the Trainer and train the model
    P = Predictor(predict_dataset=dataset,model=args.model,
        output=args.output,verbose=(not args.quiet))
    P.predict()

