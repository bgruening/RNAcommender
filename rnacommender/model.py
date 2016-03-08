import numpy as np
from theano import function, shared, config
import theano.tensor as T
from lasagne.updates import sgd

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

class Model():
    """Factorization model"""
    def __init__(self,n,m,sp,sr,irange=0.01,learning_rate=0.01,lambda_reg=0.01,seed=1234):
        """
        Params
        ------
        n : int
            Number of protein features.

        m : int
            Number of RNA features.

        sp : int
            Size of the protein latent space.

        sr : int
            Size of the RNA latent space.

        irange : float (default : 0.01)
            Initialization range for the model weights.

        learning_rate : float (default : 0.01)
            Learning rate for the weights update.

        lambda_reg : (default : 0.01)
            Lambda parameter for the regularization.

        seed : int (default : 1234)
            Seed for random number generator.
        """
        np.random.seed(seed)
        # explictit features for proteins
        fp = T.matrix("fp",dtype=config.floatX)
        # explictit features for RNAs
        fr = T.matrix("fr",dtype=config.floatX)
        # Correct label
        y = T.vector("y")

        # projection matrix for proteins
        self.ap = shared(((.5 - np.random.rand(sp,n)) * irange).astype(config.floatX), name="ap")
        self.bp = shared(((.5 - np.random.rand(sp)) * irange).astype(config.floatX), name="bp")
        # projection matrix for RNAs
        self.ar = shared(((.5 - np.random.rand(sr,m)) * irange).astype(config.floatX), name="ar")
        self.br = shared(((.5 - np.random.rand(sr)) * irange).astype(config.floatX), name="br")
        # generalization matrix
        self.B = shared(((.5 - np.random.rand(sp,sr)) * irange).astype(config.floatX), name="B")

        # Latent space for proteins
        P = T.nnet.sigmoid(T.dot(fp,self.ap.T) + self.bl)
        # Latent space for RNAs
        R = T.nnet.sigmoid(T.dot(fr,self.ar.T) + self.br)
        # Predicted output
        y_hat = T.nnet.sigmoid(T.sum(T.dot(P,self.B) * R, axis=1))

        def _regularization():
            # frobenius norm of the parameters, normalized by the size of the matrices
            norm_proteins = self.ap.norm(2) + self.bp.norm(2)
            norm_rnas = self.ar.norm(2) + self.ar.norm(2)
            norm_B = self.B.norm(2)

            num_proteins = self.ap.flatten().shape[0] + self.bp.shape[0]
            num_rnas = self.ar.flatten().shape[0] + self.br.shape[0]
            num_B = self.B.flatten().shape[0]

            return (norm_proteins/num_proteins + norm_rnas/num_rnas + norm_B/num_B)/3

        # mean squadred error
        cost_ = (T.sqr(y - y_hat)).mean()
        reg = lambda_reg*_regularization()
        cost = cost_+reg

        # compute gradient
        updates = sgd(cost, [self.w_l,self.b_l,self.w_r,self.b_r,self.C],learning_rate=learning_rate)

        # training step
        self.train = function(
                inputs=[feats_l,feats_r,y],
                outputs=[y_hat,cost],
                updates=updates)

        # predict
        self.predict = function(
                inputs=[feats_l,feats_r],
                outputs=[y_hat])

    def get_params(self):
        """Return the parameters of the model"""
        return {'ap':self.ap.get_value(),'bp':self.bp.get_value(), 'ar':self.ar.get_value(),'br':self.br.get_value(),'B':self.B.get_value()}
