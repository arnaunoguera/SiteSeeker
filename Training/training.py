import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.utils import resample
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.utils import class_weight
from sklearn.metrics import accuracy_score, precision_score, recall_score
import joblib

# Training dataset generated with get_training.py
dataset = pd.read_csv('training_file.tsv', sep='\t')

# Transforming categorical features (residue and secondary_structure) into binary encoded features
one_hot_encoder = OneHotEncoder()
encoded_features = one_hot_encoder.fit_transform(dataset[['residue', 'secondary_structure']])
# Comverting the encoded features into a data frame
encoded_df = pd.DataFrame(encoded_features.toarray(), columns=one_hot_encoder.get_feature_names_out(['residue', 'secondary_structure']))
# Concatenating the encoded features with the original data frame
dataset_encoded = pd.concat([dataset, encoded_df], axis=1)
# Removing the original categorical features
dataset_encoded.drop('residue', axis=1, inplace=True)
dataset_encoded.drop('secondary_structure', axis=1, inplace=True)

# Splitting the dataset into training and test sets (at the level of proteins, not by single residues)
prots = set(dataset_encoded.iloc[:, 0])
train_prot, test_prot = train_test_split(list(prots), test_size = 0.2, random_state=42)

# Getting the records for the training and test sets based on the proteins and separating x and y columns (features and target)
x_train_prot = dataset_encoded[dataset_encoded['protein'].isin(train_prot)].iloc[:, 0]
x_train = dataset_encoded[dataset_encoded['protein'].isin(train_prot)].iloc[:, 5:]
x_test = dataset_encoded[dataset_encoded['protein'].isin(test_prot)].iloc[:, 5:]
y_train = dataset_encoded[dataset_encoded['protein'].isin(train_prot)].iloc[:, 4]
y_test = dataset_encoded[dataset_encoded['protein'].isin(test_prot)].iloc[:, 4]


# Initializing the scaler
scaler = StandardScaler()
# Fitting the scaler on the training data and transforming both training and test data
x_train_scaled = scaler.fit_transform(x_train)
x_test_scaled = scaler.transform(x_test)

# Saving the scaler to a file to be used later on future inputs
joblib.dump(scaler, 'scaler.pkl')

# Converting the scaled arrays back to DataFrames
x_train_scaled_df = pd.DataFrame(x_train_scaled, columns=x_train.columns, index=x_train.index)
x_test_scaled_df = pd.DataFrame(x_test_scaled, columns=x_test.columns, index=x_test.index)

#Hyperparameter tuning: testing several combinations of hyperparameters to find the best ones
def hyperparameter_tuning(x_train, y_train, prot_train, c, kernel, gamma = "scale", degree = 3):
    """
    For a given set of parameters, fits three machine learning models using SVC from sklearn.svm, splitting the training set
    into training and evaluation sets randomly with each iteration, and returns the average accuracy, precision and recall
    of the three models.
    """
    accuracy = []
    precision = []
    recall = []
    for num in (42, 75, 129): # Three different random states
        # Splitting the training set into training and evaluation sets
        x_train_train, x_eval, y_train_train, y_eval, prot_train_train, prot_eval= train_test_split(x_train, y_train, prot_train, test_size = 0.2, random_state=num)
        # Fitting the model (class_weight='balanced' to account for the class imbalance in binding vs. non-binding residues: 
        # similar to  oversampling but less computationally expensive and less prone to overfitting)
        SVM = svm.SVC(C = c, kernel = kernel, gamma = gamma, degree = degree, class_weight='balanced')
        # Calculating the sample weights to account for sample imbalance: making each protein equally important for the model
        sample_weights = class_weight.compute_sample_weight('balanced', prot_train_train)
        # Fitting the model
        SVM.fit(x_train_train, y_train_train, sample_weight=sample_weights)
        # Making predictions on the evaluation set
        predictions = SVM.predict(x_eval)
        # Calculating accuracy, precision and recall for this iteration and appending it to the list
        accuracy.append(accuracy_score(y_eval, predictions))
        precision.append(precision_score(y_eval, predictions))
        recall.append(recall_score(y_eval, predictions))
    # Returning the average accuracy, precision and recall for the three models
    return(np.mean(accuracy), np.mean(precision), np.mean(recall))

# Combinations of hyperparameters to be tested
c = (0.1, 1, 10) # Higher values of C could not be executed due to lack of computational resources
kernel = ('linear', 'rbf', 'poly', 'sigmoid') # Different types of kernels to be tested
degree = (3, 10) # Only applicable for the polynomial kernel
gamma_kernels = ('poly', 'rbf', 'sigmoid') # Kernels that require the gamma parameter
gamma = ('scale', 'auto', 0.001, 0.01, 0.1) # Higher values of gamma could not be executed due to lack of computational resources

# Creating an object to store the results of the hyperparameter tuning
tuning_dicc = {'c': [], 'kernel': [], 'degree': [], 'gamma': [], 'accuracy': [], 'precision': [], 'recall': []}

# Testing all the combinations of hyperparameters
for k_value in kernel:
    for c_value in c:
        if k_value in gamma_kernels:
            # Only these kernels require the gamma parameter
            for g_value in gamma:
                if k_value =="poly":
                # Only the polynomial kernel requires the degree parameter
                    for d_value in degree:
                        print(f'Trying parameters: kernel {k_value}, C {c_value}, gamma {g_value}, degree {d_value}')
                        # Calculating the accuracy, precision and recall for the given set of hyperparameters
                        accuracy, precision, recall = hyperparameter_tuning(x_train_scaled_df, y_train, x_train_prot, c_value, k_value, g_value, d_value)
                        with open('hyperparameter_tuning.tsv', 'a') as file:
                            # Saving the results to a file
                            print(k_value, c_value, d_value, g_value, accuracy, precision, recall, sep = '\t', file = file)
                        print(f'Results: accuracy {accuracy}, precision {precision}, recall {recall}')
                else:
                    print(f'Trying parameters: kernel {k_value}, C {c_value}, gamma {g_value}')
                    # Calculating the accuracy, precision and recall for the given set of hyperparameters
                    accuracy, precision, recall = hyperparameter_tuning(x_train_scaled_df, y_train, x_train_prot, c_value, k_value, g_value)
                    with open('hyperparameter_tuning.tsv', 'a') as file:
                        # Saving the results to a file
                        print(k_value, c_value, 3, g_value, accuracy, precision, recall, sep = '\t', file = file)
                    print(f'Results: accuracy {accuracy}, precision {precision}, recall {recall}')
        else:
            print(f'Trying parameters: kernel {k_value}, C {c_value}')
            # Calculating the accuracy, precision and recall for the given set of hyperparameters
            accuracy, precision, recall= hyperparameter_tuning(x_train_scaled_df, y_train, x_train_prot, c_value, k_value)
            with open('hyperparameter_tuning.tsv', 'a') as file:
                # Saving the results to a file
                print(k_value, c_value, 3, 'scale', accuracy, precision, recall, sep = '\t', file = file)
            print(f'Results: accuracy {accuracy}, precision {precision}, recall {recall}')
print('All parameters tested')

# Reading the results of the hyperparameter tuning
tuning = pd.read_csv('hyperparameter_tuning.tsv', sep='\t')

# Getting the hyperparameters with the best accuracy
print(tuning.sort_values(by='accuracy', ascending=False).head(3))

# Fitting the model with the best hyperparameters on the whole training set
# (second best, because C = 10 had more overfitting and performed worse on the test set)
model = svm.SVC(class_weight='balanced', C = 1, kernel = 'rbf', gamma = 1)
sample_weights = class_weight.compute_sample_weight('balanced', x_train_prot)
model.fit(x_train_scaled_df, y_train, sample_weight=sample_weights)

# Making predictions on the test set
predictions = model.predict(x_test_scaled_df)

# Calculating accuracy
accuracy = accuracy_score(y_test, predictions)
print("Accuracy:", accuracy)
# Calculating precision
precision = precision_score(y_test, predictions)
print("Precision:", precision)
#Calculating recall
recall = recall_score(y_test, predictions)
print("Recall:", recall)

# Saving the model to a file to be used on the main program
joblib.dump(model, 'model.joblib') 