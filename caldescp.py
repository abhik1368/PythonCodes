import sys, cPickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn import preprocessing
 
min_max_scaler = preprocessing.MinMaxScaler()
 
trainset = sys.argv[1]
testset = sys.argv[2]
trainset = [mol for mol in Chem.SDMolSupplier(trainset) if mol is not None]
testset = [mol for mol in Chem.SDMolSupplier(testset) if mol is not None]
 
nms=[x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
print(len(nms))
 
trainDescrs = [calc.CalcDescriptors(x) for x in trainset]
testDescrs  = [calc.CalcDescriptors(x) for x in testset]
trainDescrs = np.array(trainDescrs)
testDescrs = np.array(testDescrs)
 
x_train_minmax = min_max_scaler.fit_transform( trainDescrs )
x_test_minmax = min_max_scaler.fit_transform( testDescrs )
 
classes={'(A) low':0,'(B) medium':1,'(C) high':1}
train_acts = np.array([classes[mol.GetProp("SOL_classification")] for mol in trainset],dtype="int")
test_acts = np.array([classes[mol.GetProp("SOL_classification")] for mol in testset],dtype="int")
 
dataset = ( (x_train_minmax, train_acts),(x_train_minmax, train_acts), (x_test_minmax, test_acts) )
 
f = open("rdk_sol_set_norm_descs.pkl", "wb")
cPickle.dump(dataset,f)
f.close()
