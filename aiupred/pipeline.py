from numpy import mean

import aiupred_lib
# Load the models and let AIUPred find if a GPU is available.
embedding_model, regression_model, device = aiupred_lib.init_models()
# Predict disorder of a sequence
sequence = ''
prediction = aiupred_lib.predict_disorder(sequence, embedding_model, regression_model, device)
print(prediction)
