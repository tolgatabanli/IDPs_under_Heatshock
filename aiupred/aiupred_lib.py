import logging
from scipy.signal import savgol_filter
import torch
from torch import nn, Tensor
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from torch.nn.functional import pad
import math
import os
import numpy as np

# TOLGA ADDED weights_only = True AND batch_first=True


PATH = os.path.dirname(os.path.realpath(__file__))
AA_CODE = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
WINDOW = 100


class PositionalEncoding(nn.Module):
    """
    Positional encoding for the Transformer network
    """
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        return x + self.pe[:, :x.size(1), :]


class TransformerModel(nn.Module):
    """
    Transformer model to estimate positional contact potential from an amino acid sequence
    """
    def __init__(self):
        super().__init__()
        self.d_model = 32
        self.pos_encoder = PositionalEncoding(self.d_model)
        encoder_layers = TransformerEncoderLayer(self.d_model, 2, 256, 0, batch_first=True)
        self.transformer_encoder = TransformerEncoder(encoder_layers, 2)
        self.encoder = nn.Embedding(21, self.d_model)
        self.decoder = nn.Linear((WINDOW + 1) * self.d_model, 1)

    def forward(self, src: Tensor, embed_only=False) -> Tensor:
        src = self.encoder(src) * math.sqrt(self.d_model)
        src = self.pos_encoder(src)  # (Batch x Window+1 x Embed_dim)
        embedding = self.transformer_encoder(src)
        if embed_only:
            return embedding
        output = torch.flatten(embedding, 1)
        output = self.decoder(output)
        return torch.squeeze(output)


class DecoderModel(nn.Module):
    """
    Regression model to estimate disorder propensity from and energy tensor
    """

    def __init__(self):
        super().__init__()
        input_dim = WINDOW + 1
        output_dim = 1
        current_dim = input_dim
        layer_architecture = [16, 8, 4]
        self.layers = nn.ModuleList()
        for hdim in layer_architecture:
            self.layers.append(nn.Linear(current_dim, hdim))
            current_dim = hdim
        self.layers.append(nn.Linear(current_dim, output_dim))

    def forward(self, x: Tensor) -> Tensor:
        for layer in self.layers[:-1]:
            x = torch.relu(layer(x))
        output = torch.sigmoid(self.layers[-1](x))
        return torch.squeeze(output)


@torch.no_grad()
def tokenize(sequence, device):
    """
    Tokenize an amino acid sequence. Non-standard amino acids are treated as X
    :param sequence: Amino acid sequence in string
    :param device: Device to run on. CUDA{x} or CPU
    :return: Tokenized tensors
    """
    return torch.tensor([AA_CODE.index(aa) if aa in AA_CODE else 20 for aa in sequence], device=device)


def predict_disorder(sequence, energy_model, regression_model, device, no_smoothing=False):
    """
    Predict disorder propensity from a sequence using a transformer and a regression model
    :param sequence: Amino acid sequence in string
    :param energy_model: Transformer model
    :param regression_model: regression model
    :param device: Device to run on. CUDA{x} or CPU
    :param smoothing: Use the SavGol filter to smooth the output
    :return:
    """
    predicted_energies = calculate_energy(sequence, energy_model, device)
    padded_energies = pad(predicted_energies, (WINDOW // 2, WINDOW // 2), 'constant', 0)
    unfolded_energies = padded_energies.unfold(0, WINDOW + 1, 1)
    predicted_disorder = regression_model(unfolded_energies).detach().cpu().numpy()
    if not no_smoothing and len(sequence) >= 10:
        predicted_disorder = savgol_filter(predicted_disorder, 11, 5)
    return predicted_disorder


def calculate_energy(sequence, energy_model, device):
    """
    Calculates residue energy from a sequence using a transformer network
    :param sequence: Amino acid sequence in string
    :param energy_model: Transformer model
    :param device: Device to run on. CUDA{x} or CPU
    :return: Tensor of energy values
    """
    tokenized_sequence = tokenize(sequence, device)
    padded_token = pad(tokenized_sequence, (WINDOW // 2, WINDOW // 2), 'constant', 20)
    unfolded_tokens = padded_token.unfold(0, WINDOW + 1, 1)
    return energy_model(unfolded_tokens)


def multifasta_reader(file_handler):
    """
    (multi) FASTA reader function
    :return: Dictionary with header -> sequence mapping from the file
    """
    sequence_dct = {}
    header = None
    for line in file_handler:
        if line.startswith('>'):
            header = line.strip()
            sequence_dct[header] = ''
        elif line.strip():
            sequence_dct[header] += line.strip()
    file_handler.close()
    return sequence_dct


def init_models(force_cpu=False, gpu_num=0):
    """
    Initialize networks and device to run on
    :param force_cpu: Force the method to run on CPU only mode
    :param gpu_num: Index of the GPU to use, default=0
    :return: Tuple of (embedding_model, regression_model, device)
    """
    device = torch.device(f'cuda:{gpu_num}' if torch.cuda.is_available() else 'cpu')
    if force_cpu:
        device = 'cpu'
    logging.debug(f'Running on {device}')
    if device == 'cpu':
        print('# Warning: No GPU found, running on CPU. It is advised to run AIUPred on a GPU')

    embedding_model = TransformerModel()
    embedding_model.load_state_dict(torch.load(f'{PATH}/data/embedding.pt', map_location=device, weights_only=True))
    embedding_model.to(device)
    embedding_model.eval()

    reg_model = DecoderModel()
    reg_model.load_state_dict(torch.load(f'{PATH}/data/decoder.pt', map_location=device, weights_only=True))
    reg_model.to(device)
    reg_model.eval()

    logging.debug("Networks initialized")

    return embedding_model, reg_model, device


def low_memory_predict(sequence, embedding_model, decoder_model, device, no_smoothing=False, chunk_len=1000):
    overlap = 100
    if (len(sequence)-1) % (chunk_len-overlap) == 0:
        logging.warning('Chunk length decreased by 1 to fit sequence length')
        chunk_len -= 1
    if chunk_len <= overlap:
        raise ValueError("Chunk len must be bigger than 200!")
    overlapping_predictions = []
    for chunk in range(0, len(sequence), chunk_len-overlap):
        overlapping_predictions.append(predict_disorder(
            sequence[chunk:chunk+chunk_len],
            embedding_model,
            decoder_model,
            device,
            no_smoothing
        ))
    prediction = np.concatenate((overlapping_predictions[0], *[x[overlap:] for x in overlapping_predictions[1:]]))
    return prediction


def aiupred_disorder(sequence, force_cpu=False, gpu_num=0):
    """
    Library function to carry out single sequence analysis
    :param sequence: Amino acid sequence in a string
    :param force_cpu: Force the method to run on CPU only mode
    :param gpu_num: Index of the GPU to use, default=0
    :return: Numpy array with disorder propensities for each position
    """
    embedding_model, reg_model, device = init_models(force_cpu, gpu_num)
    return predict_disorder(sequence, embedding_model, reg_model, device)


def main(multifasta_file, force_cpu=False, gpu_num=0, no_smoothing=False, low_memory=None):
    """
    Main function to be called from aiupred.py
    :param multifasta_file: Location of (multi) FASTA formatted sequences
    :param force_cpu: Force the method to run on CPU only mode
    :param gpu_num: Index of the GPU to use, default=0
    :return: Dictionary with parsed sequences and predicted results
    """
    embedding_model, reg_model, device = init_models(force_cpu, gpu_num)
    sequences = multifasta_reader(multifasta_file)
    logging.debug("Sequences read")
    logging.info(f'{len(sequences)} sequences read')
    if not sequences:
        raise ValueError("FASTA file is empty")
    results = {}
    logging.StreamHandler.terminator = ""
    for num, (ident, sequence) in enumerate(sequences.items()):
        if len(sequence) <= 10:
            logging.warning(f'Sequence length of {ident} is smaller than 10, smoothing will be turned off!\n')
        results[ident] = {}
        if low_memory:
            results[ident]['aiupred'] = low_memory_predict(sequence, embedding_model, reg_model, device, no_smoothing, low_memory)
        else:
            results[ident]['aiupred'] = predict_disorder(sequence, embedding_model, reg_model, device, no_smoothing)
        results[ident]['sequence'] = sequence
        logging.debug(f'{num}/{len(sequences)} sequences done...\r')

    logging.StreamHandler.terminator = '\n'
    logging.debug(f'Analysis done, writing output')
    return results