import os
import re
import csv
import torch
import datetime
import argparse
import pandas as pd
from tqdm import tqdm
import torch.nn.functional as F
from torch.nn import Linear
from torch.utils.data import ConcatDataset
from torch_geometric.loader import DenseDataLoader
from torch_geometric.nn import DenseGraphConv, dense_mincut_pool
from torch_geometric.data import InMemoryDataset
import torch_geometric.transforms as T


def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Train a model using TCN on spatial omics data.")
    parser.add_argument("--stage2_output_folder", type=str,
                        default="./stage2_output/", help="Output folder from stage 2.")
    parser.add_argument("--stage1_output_folder", type=str,
                        default="./stage1_output/", help="Output folder from stage 1.")
    parser.add_argument("--output_folder", type=str,
                        default="./stage3_output/", help="Output folder for this stage.")
    parser.add_argument("--image_name", type=str,
                        default="Tongue_12_10_62_TACIT_TACIT", help="Name of the image dataset.")
    parser.add_argument("--num_tcn", type=int, default=10,
                        help="Number of clusters in TCN.")
    parser.add_argument("--num_epoch", type=int, default=200,
                        help="Number of epochs for training.")
    parser.add_argument("--embedding_dimension", type=int,
                        default=128, help="Dimension of the embedding.")
    parser.add_argument("--learning_rate", type=float,
                        default=0.003, help="Learning rate for the optimizer.")
    parser.add_argument("--improvement_threshold", type=float, default=0.05,
                        help="Minimum improvement needed in six epochs to continue training.")
    parser.add_argument("--loss_cutoff", type=float, default=-0.6,
                        help="Empirical cutoff of the final loss to avoid underfitting.")
    parser.add_argument("--load_weights", action="store_true",
                        help="Load previous weights if available; otherwise, train from scratch.")
    parser.add_argument("--knn_k", type=int, default=69,
                        help="Number of neighbors to consider in KNN.")
    parser.add_argument("--save_epochs", type=int, nargs='+', default=[],
                        help="List of epochs to save the model and not delete.")
    return parser.parse_args()


class SpatialOmicsImageDataset(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super(SpatialOmicsImageDataset, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return []

    @property
    def processed_file_names(self):
        return ['SpatialOmicsImageDataset.pt']

    def download(self):
        pass


def create_datasets(stage1_output_folder, image_name, last_stage_output_folder, knn_k):
    """
    Create datasets by loading processed tiles from the previous stages.

    Args:
        stage1_output_folder (str): Path to the stage 1 output folder.
        image_name (str): Name of the image dataset.
        last_stage_output_folder (str): Path to the last stage output folder.
        knn_k (int): Number of neighbors in KNN.

    Returns:
        list: List of SpatialOmicsImageDataset objects.
    """
    datasets_paths = []
    for partition_type in ["Square", "Vertical", "Horizontal"]:
        tile_names_path = os.path.join(
            stage1_output_folder, image_name, partition_type, "tile_name_list.txt")
        tile_name_list = pd.read_csv(
            tile_names_path, sep="\t", header=None, names=["tile_name"])

        for tile_name in tile_name_list["tile_name"]:
            path = os.path.join(
                last_stage_output_folder, f"k-{knn_k}", image_name, partition_type, f"tile_{tile_name}")
            if os.path.exists(os.path.join(path, "processed")):
                datasets_paths.append(path)
            else:
                print(f"Missing processed data for {path}")

    return [SpatialOmicsImageDataset(path, transform=T.ToDense()) for path in datasets_paths]


class Net(torch.nn.Module):
    def __init__(self, in_channels, out_channels, hidden_channels=128, num_tcn=10):
        super(Net, self).__init__()
        self.conv1 = DenseGraphConv(in_channels, hidden_channels)
        self.pool1 = Linear(hidden_channels, num_tcn)

    def forward(self, x, adj, mask=None):
        x = F.relu(self.conv1(x, adj, mask))
        s = self.pool1(x)
        x, adj, mc1, o1 = dense_mincut_pool(x, adj, s, mask)
        return F.log_softmax(x, dim=-1), mc1, o1, s, adj


def train_epoch(model, train_loader, optimizer, device, accumulation_steps=32):
    """
    Train the model for one epoch.

    Args:
        model (torch.nn.Module): The model to train.
        train_loader (DataLoader): DataLoader for the training data.
        optimizer (torch.optim.Optimizer): The optimizer to use.
        device (torch.device): The device to run on (CPU or GPU).
        accumulation_steps (int): Number of steps to accumulate gradients before updating.

    Returns:
        float: Average loss over the epoch.
    """
    model.train()
    loss_all = 0
    optimizer.zero_grad()

    for i, data in enumerate(train_loader):
        data = data.to(device)
        out, mc_loss, o_loss, _, _ = model(data.x, data.adj, data.mask)
        loss = (mc_loss + o_loss) / accumulation_steps
        loss.backward()

        if (i + 1) % accumulation_steps == 0:
            optimizer.step()
            optimizer.zero_grad()

        loss_all += loss.item() * accumulation_steps

    if (i + 1) % accumulation_steps != 0:
        optimizer.step()
        optimizer.zero_grad()

    return loss_all / len(train_loader)


def get_latest_checkpoint(run_folder):
    """
    Get the latest checkpoint from the run folder.

    Args:
        run_folder (str): Path to the run folder.

    Returns:
        tuple: Path to the latest checkpoint and the epoch number.
    """
    checkpoints = [f for f in os.listdir(
        run_folder) if re.match(r'TCN_model_epoch\d+\.pt', f)]
    if not checkpoints:
        return None, 0

    def extract_epoch(filename):
        match = re.search(r'epoch(\d+)', filename)
        return int(match.group(1)) if match else -1

    latest_checkpoint = max(checkpoints, key=extract_epoch)
    epoch_number = extract_epoch(latest_checkpoint)

    return os.path.join(run_folder, latest_checkpoint), epoch_number


def save_model_at_epoch(model, epoch, run_folder_name, saved_models, save_epochs):
    """
    Save the model at a given epoch and manage the saved models.

    Args:
        model (torch.nn.Module): The model to save.
        epoch (int): The current epoch number.
        run_folder_name (str): The folder where the model is saved.
        saved_models (list): List of saved model filenames.
        save_epochs (list): List of epochs at which the model should not be deleted.
    """
    model_filename = os.path.join(
        run_folder_name, f"TCN_model_epoch{epoch}.pt")
    torch.save(model.state_dict(), model_filename)

    if epoch not in save_epochs:
        saved_models.append(model_filename)
        if len(saved_models) > 3:
            oldest_model = saved_models.pop(0)
            if os.path.exists(oldest_model):
                os.remove(oldest_model)


def load_weights_if_available(args, model, run_folder_name):
    """
    Load weights from the latest checkpoint if available.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        model (torch.nn.Module): The model to load weights into.
        run_folder_name (str): The folder where the model is saved.

    Returns:
        int: The starting epoch number.
    """
    start_epoch = 1
    if args.load_weights:
        latest_checkpoint, start_epoch = get_latest_checkpoint(run_folder_name)
        if latest_checkpoint:
            print(
                f"Loading weights from {latest_checkpoint} starting at epoch {start_epoch + 1}")
            model.load_state_dict(torch.load(latest_checkpoint))
            start_epoch += 1
        else:
            print("No previous weights found. Starting training from scratch.")
    return start_epoch


def initialize_loss_file(filename, start_epoch, load_weights):
    """
    Initialize the loss CSV file.

    Args:
        filename (str): The filename for the loss CSV.
        start_epoch (int): The starting epoch number.
        load_weights (bool): Whether to load weights from a previous run.

    Returns:
        None
    """
    headers = ["Epoch", "UnsupervisedLoss"]
    file_mode = 'a' if load_weights and start_epoch > 1 else 'w'
    with open(filename, file_mode, newline='') as f:
        f_csv = csv.writer(f)
        if file_mode == 'w':
            f_csv.writerow(headers)


def log_epoch_loss(filename, epoch, loss):
    """
    Log the loss for an epoch.

    Args:
        filename (str): The filename for the loss CSV.
        epoch (int): The epoch number.
        loss (float): The loss value.

    Returns:
        None
    """
    with open(filename, "a", newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow([epoch, loss])


def should_stop_early(epoch, train_loss, previous_losses, args):
    """
    Determine if training should stop early due to lack of improvement.

    Args:
        epoch (int): The current epoch number.
        train_loss (float): The training loss for the current epoch.
        previous_losses (list): List of previous losses.
        args (argparse.Namespace): Parsed command-line arguments.

    Returns:
        bool: True if training should stop early, False otherwise.
    """
    return all(abs(prev_loss - train_loss) < args.improvement_threshold for prev_loss in previous_losses) \
        and epoch > args.num_epoch / 2 and min(previous_losses) < args.loss_cutoff


def update_previous_losses(previous_losses, train_loss):
    """
    Update the list of previous losses.

    Args:
        previous_losses (list): List of previous losses.
        train_loss (float): The training loss for the current epoch.

    Returns:
        None
    """
    previous_losses.pop(0)
    previous_losses.append(train_loss)


def train_and_evaluate_model(args, datasets, output_folder_name):
    """
    Train and evaluate the model.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        datasets (list): List of datasets to train on.
        output_folder_name (str): The folder to save the outputs.

    Returns:
        None
    """
    train_dataset = ConcatDataset(datasets)
    train_loader = DenseDataLoader(train_dataset, batch_size=1, shuffle=True)

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = Net(datasets[0].num_features, 1,
                hidden_channels=args.embedding_dimension, num_tcn=args.num_tcn).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)

    run_folder_name = output_folder_name
    os.makedirs(run_folder_name, exist_ok=True)

    start_epoch = load_weights_if_available(args, model, run_folder_name)

    filename_0 = os.path.join(run_folder_name, "unsupervised_epochs.csv")
    initialize_loss_file(filename_0, start_epoch, args.load_weights)

    previous_losses = [float("inf")] * 6
    saved_models = []

    with tqdm(total=args.num_epoch, initial=start_epoch-1, desc="Training Progress", unit="epoch") as pbar:
        for epoch in range(start_epoch, args.num_epoch + 1):
            train_loss = train_epoch(model, train_loader, optimizer, device)
            pbar.set_postfix({"Train Loss": train_loss})
            pbar.update(1)

            log_epoch_loss(filename_0, epoch, train_loss)

            if should_stop_early(epoch, train_loss, previous_losses, args):
                save_model_at_epoch(
                    model, epoch, run_folder_name, saved_models, args.save_epochs)
                break

            if epoch % 20 == 0 or epoch in args.save_epochs:
                save_model_at_epoch(
                    model, epoch, run_folder_name, saved_models, args.save_epochs)

            update_previous_losses(previous_losses, train_loss)

    print(f"Final train loss is {train_loss:.4f}")


def main(args):
    datasets = create_datasets(
        args.stage1_output_folder, args.image_name, args.stage2_output_folder, args.knn_k)
    output_folder_name = os.path.join(
        args.output_folder, f"k-{args.knn_k}_TCN-{args.num_tcn}_lr-{args.learning_rate}", args.image_name)
    os.makedirs(output_folder_name, exist_ok=True)

    train_and_evaluate_model(args, datasets, output_folder_name)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
