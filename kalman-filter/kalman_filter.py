#!/usr/bin/python

# Import python libraries
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt


def Rz(q):
    """Rotation matrix around z"""
    return np.array([[np.cos(q), -np.sin(q), 0],
                     [np.sin(q), np.cos(q), 0],
                     [0, 0, 1]])

def plot_ellipse(ax, center, cov_matrix):
    w, R = np.linalg.eig(cov_matrix)
    S = np.diag(np.sqrt(w))
    Transf = R.dot(S)

    q = np.linspace(0, 2*np.pi, 100)
    for th in np.linspace(-np.pi/2, np.pi/2, 8):
        # horizontal lines
        radius = np.cos(th)
        height = np.sin(th)
        xyz = np.stack((radius*np.cos(q), radius*np.sin(q), height*np.ones(q.shape))).T
        pts = xyz.dot(Transf.T) + center.reshape(3)
        plt.plot(pts[:,0], pts[:,1], pts[:,2], 'k', lw=0.5, alpha=0.3)

        # vertical lines
        xyz = Rz(th).dot(np.stack((np.cos(q), np.zeros(q.shape), np.sin(q)))).T
        pts = xyz.dot(Transf.T) + center.reshape(3)
        plt.plot(pts[:,0], pts[:,1], pts[:,2], 'k', lw=0.5, alpha=0.3)



class KalmanFilter(object):
    """Kalman Filter class keeps track of the estimated state of
    the system and the variance or uncertainty of the estimate.
    Predict and Correct methods implement the functionality
    Reference: https://en.wikipedia.org/wiki/Kalman_filter
    Attributes: None
    """

    def __init__(self, F=1.0, B=1.0, H=1.0, Q=1.0, R=1.0, X=0.0, P=1.0):
        """Initialize variable used by Kalman Filter class
        Args:
            F - state transition matrix
            B - contron input matrix
            H - measurement matrix
            Q - process noise covariance
            R - measurement noise covariance
            X - initial state
            P - state error covariance
        """
        (self.F, self.B, self.H, self.Q, self.R, self.X, self.P) = (F, B, H, Q, R, X, P)

    def predict(self, U):
        """Predict state vector X and variance of uncertainty P (covariance).
        Args:
            U - control input
        Return:
            vector of predicted state estimate and covariance
        """
        # Predicted state estimate
        self.X = self.F.dot(self.X) + self.B.dot(U)
        self.P = self.F.dot(self.P.dot(self.F.T)) + self.Q

        return self.X, self.P

    def correct(self, Z):
        """Correct or update state vector u and variance of uncertainty P (covariance).
        Args:
            Z: vector of observations
        Return:
            vector of predicted state estimate and covariance
        """

        Yt = Z - self.H.dot(self.X)
        S = self.H.dot(self.P.dot(self.H.T)) + self.R  # cov of Y priori
        K = np.linalg.lstsq(S.T, self.H.dot(self.P))[0].T
        self.X += K.dot(Yt)
        self.P = (np.eye(self.X.shape[0]) - K.dot(self.H)).dot(self.P)
        Yt = Z - self.H.dot(self.X)

        return self.X, self.P


class SimpleUAV3D(KalmanFilter):
    """
    TODO
     - For new detections, initial velocity guess should be towards the camera (relatively)
     - [Validation gate for measurements](http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf)
     
    """ 

    def __init__(self, m=1.0, T=1.0, X=np.zeros((6,1)), P=np.eye(6), Q=np.eye(6), R=np.eye(3)):

        F = np.array([[1, 0, 0, T, 0, 0],
                    [0, 1, 0, 0, T, 0],
                    [0, 0, 1, 0, 0, T],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1]], np.float)  # matrix in observation equations
        B = np.array([[T/2, 0, 0],
                    [0, T/2, 0],
                    [0, 0, T/2],
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]], np.float) * T/m  # previous state vector
        H = np.array([[1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0]], np.float)  # matrix in observation equations
        super(SimpleUAV3D, self).__init__(F, B, H, Q, R, X, P)

    def plot_state(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')

        ax.cla()
        position = self.X[:3]
        velocity = self.X[3:6]
        ax.quiver(position[0], position[1], position[2], 
                    velocity[0], velocity[1], velocity[2],
                    length=1, normalize=True)
        plot_ellipse(ax, position, self.P[:3,:3])

        ax.set_xlabel('X'), ax.set_ylabel('Y'), ax.set_zlabel('Z')
        ax.set_xlim(-10,10), ax.set_ylim(-10,10), ax.set_zlim(-10,10)
        ax.set_aspect('equal', 'box')

        return ax

