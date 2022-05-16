import numpy as np

from scipy.optimize import fmin_cobyla

import openTSNE


def cartesian2polar(X,Y):
    """
    Convert cartsian coordinates to polar (https://www.varsitytutors.com/hotmath/hotmath_help/topics/polar-coordinates)
    """
    r = np.sqrt(X**2 + Y**2)
    t = np.arctan2(Y, X)
    return r, t

def translate_function(X,Y, x_center = 0, y_center = 0):
    """
    Translate function origin point
    """
    
    X_trasl = X - x_center
    Y_trasl = Y - y_center
    
    return X_trasl, Y_trasl

def boundary_function_cartesian(X,Y):
    """
    Write a generic function F in Cartesian coordinates'
    The example below plots a circle of radius = 1
    """
    
    F = X**2 + Y**2 - 1 #circle
    
    return F

def boundary_circle(X, Y, radius = 1, x_center = 0, y_center = 0, **kwargs):
    
    # translate circle
    X_trasl, Y_trasl = translate_function(X,Y, x_center, y_center)
    
    # convert Cartesian 
    r , _ =  cartesian2polar(X_trasl,Y_trasl)
    
    # circle function in polar coordinates
    F = radius - r
    
    return F

def boundary_heart(X, Y, size = 1, x_center = 0, y_center = 0, **kwargs):
    
    # translate circle
    X_trasl, Y_trasl = translate_function(X,Y, x_center, y_center)
    
    # convert Cartesian 
    r , t =  cartesian2polar(X_trasl,Y_trasl)
    
    # heart function in polar coordinates (https://mathworld.wolfram.com/HeartCurve.html)
    F = np.sin(t) * np.sqrt(np.absolute(np.cos(t))) / (np.sin(t) + 7/5) - 2*np.sin(t) - r/size + 2
    
    return F

def closest_point_in_function(point, function):
    def objective(p_close):
        x,y = p_close
        return np.sqrt((x - point[0])**2 + (y - point[1])**2)

    def constrain(p_close):
        x,y = p_close
        return function(x,y)

    p = fmin_cobyla(objective, x0=point, cons=[constrain])
    
    return p

def circular_openTSNE(adata, radius = 1, x_center = 0, y_center = 0):
    
    buondary_function = lambda x , y : boundary_circle(x, y, radius = radius, x_center = x_center, y_center = y_center)

    def set_boundary_embedding(emb, embeddings, boundary_function):
        outside_boundary = buondary_function(emb[:,0], emb[:,1]) < 0
        
        if any(outside_boundary):
            output = list(map(lambda point : closest_point_in_function(point, buondary_function), emb[outside_boundary]))
            emb[outside_boundary] = np.array(output)
        
        embeddings.append(np.array(emb))
    
    if adata.shape[0]==1:
        adata.obsm['X_tsne'] = np.array([[x_center, y_center]])
    else:
        embeddings = []

        tsne = openTSNE.TSNE(perplexity=50, metric="cosine", n_jobs=32, verbose=True,
                             # The embedding will be appended to the list we defined above, make sure we copy the
                             # embedding, otherwise the same object reference will be stored for every iteration
                             callbacks=lambda it, err, emb: set_boundary_embedding(emb, embeddings, buondary_function),
                             # This should be done on every iteration
                             callbacks_every_iters=1)

        tsne_embedding = tsne.fit(adata.obsm['X_pca'])

        adata.obsm['X_tsne'] = embeddings[-1]