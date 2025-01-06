
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
# x(t+ h) = x(t) + h * ( x(t) * alpha - x(t) * beta * y(t))
# y(t+ h) = y(t) + h * ( y(t) * delta * x(t) - y(t) * gamma) 


def euler(lapin ,renard, h = 0.001, alpha = 1/3, beta =1/3, delta =1/3, gamma =1/3 ) :

    lapin_euler = [lapin[0] / 1000]
    renard_euler = [renard[0] / 1000]

    #nombre_observations - 1

    nombre_observations = len(lapin)
    for i in range(0, nombre_observations - 1):
        
        next_value_lapin = lapin_euler[i] + h * lapin_euler[i] * (alpha - beta * renard_euler[i])
        next_value_renard = renard_euler[i] + h * renard_euler[i] * (delta * lapin_euler[i] - gamma)

        

        lapin_euler.append(next_value_lapin)
        renard_euler.append(next_value_renard)
        # print('effectifs lapins : {}'.format(next_value_lapin))
        # print('effectifs renard : {}'.format(next_value_renard))"""
    return lapin_euler, renard_euler

def mse(lapin_real, renard_real, lapin_euler, renard_euler):
    lapin_real = lapin_real / 1000
    renard_real = renard_real / 1000

    square_error_lapin = []
    square_error_renard = []
    for i in range(len(lapin_real)):
        square_error_lapin.append((lapin_real[i] - lapin_euler[i])**2)
        square_error_renard.append((renard_real[i] - renard_euler[i])**2)
    mse_lapin = sum(square_error_lapin)/len(square_error_lapin)
    mse_renard = sum(square_error_renard)/len(square_error_renard)

    return mse_lapin + mse_renard



def gride_search(lapin, renard, h , param):
    param_values = param.values()
    combinations = list(itertools.product(*param_values))

    best_score = float('inf')  
    best_params = None
    best_lapins, best_renards = None, None

    for combination in combinations:
        alpha, beta, delta, gamma = combination
        lapin_euler, renard_euler = euler(data['lapin'], data['renard'], h, alpha, beta, delta, gamma)
        mse_global = mse(data['lapin'], data['renard'], lapin_euler, renard_euler)

    

        if mse_global < best_score:
            best_score = mse_global
            best_params = combination
            best_lapins, best_renards = lapin_euler, renard_euler

    return best_params, best_score, best_lapins, best_renards


if __name__ == "__main__":
    data = pd.read_csv('populations_lapins_renards.csv')
    lapin_euler, renard_euler = euler(data['lapin'], data['renard'], h=0.001, alpha= 2/3, beta=4/3, delta=1, gamma=1)
    mse_global = mse(data['lapin'], data['renard'],lapin_euler, renard_euler)
    print('mse globale : {}'.format(np.round(mse_global, 3)))

    plt.figure(figsize=(15, 6))
    plt.plot(data['date'],np.array(lapin_euler), ":")
    plt.plot(data['date'],np.array(renard_euler), color="red")
    plt.show()  

    param_dict = {
    'alpha': [1/3, 2/3, 1, 4/3],
    'beta': [1/3, 2/3, 1, 4/3],
    'delta': [1/3, 2/3, 1, 4/3],
    'gamma': [1/3, 2/3, 1, 4/3],
    } 

    best_params, best_score, best_lapins, best_renards = gride_search(data['lapin'], data['renard'], 0.001 , param_dict)
    print("Meilleurs paramètres trouvés:", best_params)
    print("Meilleur score (RMSE):", best_score) 

    plt.figure(figsize=(15, 6))
    plt.plot(data['date'],best_lapins, ":")
    plt.plot(data['date'],best_renards, color="red")
    plt.show()

    
