import numpy as np
from numpy.linalg import inv,pinv,svd

def dem_inv_gsvd(A,B):
    """
    dem_inv_gsvd

    Performs the generalised singular value decomposition of two matrices A,B.

    Inputs

    A:
        cross section matrix
    B:
        regularisation matrix (square)

    Performs

    the decomposition of:

        A=U*SA*W^-1
        B=V*SB*W^-1

        with gsvd matrices u,v and the weight W and diagnoal matrics SA and SB

    Outputs

    U:
        decomposition product matrix
    V:
        decomposition prodyct matrix
    W:
        decomposition prodyct matrix
    alpha:
        the vector of the diagonal values of SA
    beta:
        the vector of the diagonal values of SB
  

    """  
    #calculate the matrix A*B^-1
    AB1=A@inv(B)
    sze=AB1.shape
    C=np.zeros([max(sze),max(sze)])
    C[:sze[0],:sze[1]]=AB1
    #use np.linalg.svd to calculate the singular value decomposition
    u,s,v = svd(C,full_matrices=True,compute_uv=True)
    # U, S, Vh = svd(AB1, full_matrices=False)
    #from the svd products calculate the diagonal components form the gsvd
    beta=1./np.sqrt(1+s**2)
    alpha=s*beta

    #diagonalise alpha and beta into SA and SB
    #onea=np.diag(alpha)
    oneb=np.diag(beta)
     #calculate the w matrix
    # w=inv(inv(onea)@transpose(u)@A)
    w2=pinv(inv(oneb)@v@B)

    #return gsvd products, transposing v as we do.
    return alpha,beta,u.T[:,:sze[0]],v.T,w2
