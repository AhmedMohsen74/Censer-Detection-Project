import itertools
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from pandas.plotting import scatter_matrix
import seaborn as sns
import pickle

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
names = ['id', 'clump_thickness', 'uniform_cell_size', 'uniform_cell_shape',
         'marginal_adhesion', 'single_epithelial_size', 'bare_nuclei',
         'bland_chromatin', 'normal_nucleoli', 'mitoses', 'class']
df = pd.read_csv(url, names=names)
df.drop(['id'], axis=1, inplace=True)
df.replace('?', np.nan, inplace=True)
df.fillna(method='ffill', inplace=True)
df['bare_nuclei'] = df['bare_nuclei'].astype('int64')
sns.displot(df['class'], kde=True)
ax = df[df['class'] == 4][0:50].plot(
    kind='scatter', x='clump_thickness', y='uniform_cell_size', color='DarkBlue', label='malignant')
df[df['class'] == 2][0:50].plot(kind='scatter', x='clump_thickness',
                                y='uniform_cell_size', color='Yellow', label='benign', ax=ax)
plt.show()
sns.set_style('darkgrid')
df.hist(figsize=(30, 30))
plt.show()
scatter_matrix(df, figsize=(18, 18))
plt.show()
plt.figure(figsize=(10, 10))
sns.boxplot(data=df, orient='h')
df.corr()
plt.figure(figsize=(30, 20))
cor = df.corr()
sns.heatmap(cor, vmax=1, square=True, annot=True, cmap=plt.cm.Blues)
plt.title('Correlation between different attributes')
plt.show()
sns.pairplot(df, diag_kind='kde')

# Correlation with output variable
cor_target = abs(cor["class"])
# Selecting highly correlated features
relevant_features = cor_target[cor_target > 0]
# relevant_features

# Split the data into predictor variables and target variable, following by breaking them into train and test sets.

Y = df['class'].values
X = df.drop('class', axis=1).values

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.30, random_state=21)

# Testing Options
scoring = 'accuracy'

# Define models to train
models = []
models.append(('CART', DecisionTreeClassifier()))
models.append(('SVM', SVC()))
models.append(('NB', GaussianNB()))
models.append(('KNN', KNeighborsClassifier()))

# evaluate each model in turn
results = []
names = []

for name, model in models:
    kfold = KFold(n_splits=10)
    cv_results = cross_val_score(
        model, X_train, Y_train, cv=kfold, scoring=scoring)
    results.append(cv_results)
    names.append(name)
    msg = "For %s Model:Mean accuracy is %f (Std accuracy is %f)" % (
        name, cv_results.mean(), cv_results.std())
    print(msg)

fig = plt.figure(figsize=(10, 10))
fig.suptitle('Performance Comparison')
ax = fig.add_subplot(111)
plt.boxplot(results)
ax.set_xticklabels(names)
plt.show()

# Make predictions on validation dataset

for name, model in models:
    model.fit(X_train, Y_train)
    predictions = model.predict(X_test)
    print("\nModel:", name)
    print("Accuracy score:", accuracy_score(Y_test, predictions))
    print("Classification report:\n", classification_report(Y_test, predictions))


clf = SVC()

clf.fit(X_train, Y_train)
accuracy = clf.score(X_test, Y_test)
print("Test Accuracy:", accuracy)

predict = clf.predict(X_test)
# predict

example_measures = [[4, 2, 1, 1, 1, 2, 3, 2, 1]]
prediction = clf.predict(example_measures)
print(prediction)

sns.set_theme(style="dark")


def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')


# Compute confusion matrix
cnf_matrix = confusion_matrix(Y_test, predict, labels=[2, 4])
np.set_printoptions(precision=2)

print(classification_report(Y_test, predict))

# Plot non-normalized confusion matrix
plt.figure()
plot_confusion_matrix(cnf_matrix, classes=[
                      'Benign(2)', 'Malignant(4)'], normalize=False,  title='Confusion matrix')

pickle.dump(clf, open('model.pkl', 'wb'))

model = pickle.load(open('model.pkl', 'rb'))
print(model.predict([[4, 2, 1, 1, 1, 2, 3, 2, 1]]))


st.title("Detecting Breast Cancer")
st.sidebar.title('Breast Cancer Detection')
navigation = st.sidebar.radio('VIEW', ('Data Analysis', 'Prediction'))

seq1 = st.number_input("Clump Thikness: ", key="seq1")
seq2 = st.number_input("Uniform Cell Size: ", key="seq2")
seq3 = st.number_input("Uniform Cell Shape: ", key="seq3")
seq4 = st.number_input("Marginal Adhesion: ", key="seq4")
seq5 = st.number_input("Single Epithelial Cell Size: ", key="seq5")
seq6 = st.number_input("Bare Nuclei: ", key="seq6")
seq7 = st.number_input("Bland Chromatin: ", key="seq7")
seq8 = st.number_input("Normal Nucleoli: ", key="seq8")
seq9 = st.number_input("Mitosis: ", key="seq9")


if(st.button("Predict Cancer")):
    input_features = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9]
    features_value = [np.array(input_features)]
    features_name = ['clump_thickness', 'uniform_cell_size', 'uniform_cell_shape',
                     'marginal_adhesion', 'single_epithelial_size', 'bare_nuclei',
                     'bland_chromatin', 'normal_nucleoli', 'mitoses']
    df = pd.DataFrame(features_value, columns=features_name)
    output = model.predict(df)
    if output == 4:
        st.write("Patient Has Breast cancer")
    else:
        st.write("Patient Has No Breast cancer")
