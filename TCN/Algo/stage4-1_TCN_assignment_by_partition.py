import os
import re
import torch
import datetime
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch.nn.functional as F
from torch.nn import Linear
from torch.utils.data import ConcatDataset
from torch_geometric.loader import DenseDataLoader
from torch_geometric.data import InMemoryDataset
from torch_geometric.nn import DenseGraphConv, dense_mincut_pool
import torch_geometric.transforms as T


def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Predict TCN using a trained TCN model by partition methods.")
    parser.add_argument("--stage1_output_folder", type=str,
                        default="./stage1_output/", help="Output folder from stage 1.")
    parser.add_argument("--stage2_output_folder", type=str,
                        default="./stage2_output/", help="Output folder from stage 2.")
    parser.add_argument("--stage3_output_folder", type=str,
                        default="./stage3_output/", help="Output folder from stage 3.")
    parser.add_argument("--output_folder", type=str,
                        default="./stage4-1_output/", help="Output folder.")
    parser.add_argument("--load_config", type=str,
                        default="k-69_TCN-10_lr-0.003/", help="Configuration of model to load.")
    parser.add_argument("--image_name", type=str, default="Tongue_12_10_62_TACIT_TACIT",
                        help="Name of the image dataset to be processed.")
    parser.add_argument("--num_tcn", type=int, default=10,
                        help="Number of clusters in TCN.")
    parser.add_argument("--embedding_dimension", type=int,
                        default=128, help="Dimension of the embedding.")
    parser.add_argument("--knn_k", type=int, default=69,
                        help="Number of neighbors to consider in KNN.")
    parser.add_argument("--load_epoch", type=int, default=None,
                        help="Specific epoch to load the model from. If not given, the latest checkpoint is used.")
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


def create_datasets(stage1_output_folder, image_name, stage2_output_folder, knn_k):
    """
    Create datasets by loading processed tiles from the previous stages.

    Args:
        stage1_output_folder (str): Path to the stage 1 output folder.
        image_name (str): Name of the image dataset.
        stage2_output_folder (str): Path to the stage 2 output folder.
        knn_k (int): Number of neighbors in KNN.

    Returns:
        tuple: Lists of dataset paths for Square, Vertical, and Horizontal partitions.
    """
    dataset_paths = {
        "Square": [],
        "Vertical": [],
        "Horizontal": []
    }

    for partition_type in ["Square", "Vertical", "Horizontal"]:
        tile_name_list_path = os.path.join(
            stage1_output_folder, image_name, partition_type, "tile_name_list.txt")
        tile_name_list = pd.read_csv(
            tile_name_list_path, sep="\t", header=None, names=["tile_name"])

        for tile_name in tile_name_list["tile_name"]:
            path = os.path.join(
                stage2_output_folder, f"k-{knn_k}", image_name, partition_type, f"tile_{tile_name}")
            if os.path.exists(os.path.join(path, "processed")):
                dataset_paths[partition_type].append(path)
            else:
                print(f"Missing processed data for {path}")

    return dataset_paths["Square"], dataset_paths["Vertical"], dataset_paths["Horizontal"]


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


def get_checkpoint_by_epoch(folder, epoch):
    """
    Get the checkpoint file for a specific epoch.

    Args:
        folder (str): Folder containing checkpoint files.
        epoch (int): The specific epoch to load.

    Returns:
        str: Path to the checkpoint file.
    """
    checkpoint_file = f"TCN_model_epoch{epoch}.pt"
    if os.path.exists(os.path.join(folder, checkpoint_file)):
        return os.path.join(folder, checkpoint_file)
    return None


def get_latest_checkpoint(folder):
    """
    Get the latest checkpoint file from a folder.

    Args:
        folder (str): Folder containing checkpoint files.

    Returns:
        tuple: Path to the latest checkpoint file and the epoch number.
    """
    checkpoint_files = [f for f in os.listdir(
        folder) if re.match(r'TCN_model_epoch\d+\.pt', f)]
    if not checkpoint_files:
        return None, None

    latest_epoch = max(int(re.findall(r'\d+', f)[0]) for f in checkpoint_files)
    return os.path.join(folder, f"TCN_model_epoch{latest_epoch}.pt"), latest_epoch


def process_partition_data(loader, model, run_folder, partition_type, device):
    """
    Process and save the results for a partition.

    Args:
        loader (DataLoader): DataLoader for the partition data.
        model (torch.nn.Module): The model to use for processing.
        run_folder (str): Folder to save the results.
        partition_type (str): Type of partition (Square, Vertical, Horizontal).
        device (torch.device): Device to run the model on.
    """
    for data in tqdm(loader, desc=f"Processing {partition_type} Data", unit="batch"):
        data = data.to(device)
        model.eval()
        result = model(data.x, data.adj, data.mask)

        cluster_assign_matrix = torch.softmax(
            result[3][0, :, :], dim=-1).detach().cpu().numpy()
        cluster_adj_matrix = result[4][0, :, :].detach().cpu().numpy()
        node_mask = data.mask.detach().cpu().numpy().T

        with open(os.path.join(run_folder, f"{partition_type}_TCN_assign_matrix.csv"), 'a') as f:
            np.savetxt(f, cluster_assign_matrix, delimiter=',')

        with open(os.path.join(run_folder, f"{partition_type}_TCN_adjacent_matrix.csv"), 'a') as f:
            np.savetxt(f, cluster_adj_matrix, delimiter=',')

        with open(os.path.join(run_folder, f"{partition_type}_node_mask.csv"), 'a') as f:
            np.savetxt(f, node_mask, delimiter=',', fmt='%i')



def main(args):
    square_paths, vertical_paths, horizontal_paths = create_datasets(
        args.stage1_output_folder, args.image_name, args.stage2_output_folder, args.knn_k
    )

    datasets = {
        "Square": ConcatDataset([SpatialOmicsImageDataset(path, transform=T.ToDense()) for path in square_paths]),
        "Vertical": ConcatDataset([SpatialOmicsImageDataset(path, transform=T.ToDense()) for path in vertical_paths]),
        "Horizontal": ConcatDataset([SpatialOmicsImageDataset(path, transform=T.ToDense()) for path in horizontal_paths])
    }

    loaders = {
        "Square": DenseDataLoader(datasets["Square"], batch_size=1, shuffle=False),
        "Vertical": DenseDataLoader(datasets["Vertical"], batch_size=1, shuffle=False),
        "Horizontal": DenseDataLoader(datasets["Horizontal"], batch_size=1, shuffle=False)
    }

    run_folder = os.path.join(
        args.stage3_output_folder, args.load_config, args.image_name)

    if not os.path.exists(run_folder):
        print(f"Run folder {run_folder} does not exist.")
        return

    if args.load_epoch is not None:
        checkpoint = get_checkpoint_by_epoch(run_folder, args.load_epoch)
        epoch_to_use = args.load_epoch
    else:
        checkpoint, epoch_to_use = get_latest_checkpoint(run_folder)

    if not checkpoint:
        print("No checkpoint found.")
        return

    print(f"Loading model from {checkpoint}")

    output_folder = os.path.join(
        args.output_folder, args.load_config, args.image_name, f"checkpoint_epoch_{epoch_to_use}")
    os.makedirs(output_folder, exist_ok=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = Net(datasets["Square"][0].num_features, 1,
                hidden_channels=args.embedding_dimension, num_tcn=args.num_tcn)
    model.load_state_dict(torch.load(checkpoint))
    model.to(device)

    for partition_type, loader in loaders.items():
        process_partition_data(
            loader, model, output_folder, partition_type, device)

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
