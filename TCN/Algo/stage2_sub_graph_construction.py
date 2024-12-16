import argparse
import os
import shutil
import datetime
import numpy as np
import pandas as pd
import torch
from sklearn.neighbors import kneighbors_graph
import torch_geometric.transforms as T
from torch_geometric.data import Data, InMemoryDataset


def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate topology structures and node attribute matrices for KNN graphs of every tile created in Stage 1.")
    parser.add_argument("--data_source", type=str, default="./PCF_data/",
                        help="Root directory containing image data files.")
    parser.add_argument("--input_folder", type=str, default="./stage0_output/",
                        help="Input folder containing partitioned data.")
    parser.add_argument("--knn_k", type=int, default=69,
                        help="Number of neighbors to consider in KNN.")
    parser.add_argument("--output_folder", type=str,
                        default="./stage1_output/", help="Output folder to save results.")
    return parser.parse_args()


def load_region_names(data_source):
    """
    Load region names from the data source file.

    Args:
        data_source (str): Path to the data source directory.

    Returns:
        pd.DataFrame: DataFrame containing region names.
    """
    region_filename = os.path.join(data_source, "sample_name_list.txt")
    return pd.read_csv(region_filename, sep="\t", header=None, names=["slide"])


def prepare_output_folder(output_folder, knn_k, reset=False):
    """
    Prepare the output directory, removing it first if it already exists.

    Args:
        output_folder (str): Root directory for output.
        knn_k (int): Number of neighbors in KNN.

    Returns:
        str: Path to the prepared output folder.
    """
    output_folder_knn = os.path.join(output_folder, f"k-{knn_k}")
    if os.path.exists(output_folder_knn) and reset:
        shutil.rmtree(output_folder_knn)
        os.makedirs(output_folder_knn)
    return output_folder_knn


def process_each_tile(tile_name_list, tile_root, knn_k, output_folder, region_name, partition_type):
    """
    Process each tile to generate and save KNN graphs.

    Args:
        tile_name_list (pd.DataFrame): List of tile names.
        tile_root (str): Root directory of the tiles.
        knn_k (int): Number of neighbors in KNN.
        output_folder (str): Directory to save the results.
        region_name (str): Name of the region.
        partition_type (str): Type of partition (Square, Vertical, Horizontal).

    Returns:
        int: Number of points skipped in the last tile.
    """
    accumulated_coords = np.empty((0, 2))
    skipped_points_count = 0

    for tile_name in tile_name_list["tile_name"]:
        output_tile_dir = os.path.join(
            output_folder, region_name, partition_type, f"tile_{tile_name}")
        os.makedirs(output_tile_dir, exist_ok=True)

        graph_coord_filename = os.path.join(
            tile_root, f"tile_{tile_name}_coordinate.txt")
        x_y_coordinates = np.loadtxt(
            graph_coord_filename, dtype='float', delimiter="\t").reshape(-1, 2)

        if accumulated_coords.size > 0:
            x_y_coordinates = np.vstack((accumulated_coords, x_y_coordinates))
            accumulated_coords = np.empty((0, 2))

        if x_y_coordinates.shape[0] < knn_k:
            print(
                f"Sample-{region_name} tile-{tile_name} has less than {knn_k} cells. Moving to next tile!")
            accumulated_coords = x_y_coordinates if tile_name_list.index[-1] > tile_name_list.index[tile_name] else np.empty(
                (0, 2))
            skipped_points_count = x_y_coordinates.shape[0] if accumulated_coords.size == 0 else 0
            continue

        knn_graph = kneighbors_graph(
            x_y_coordinates, knn_k, mode='connectivity', include_self=False, n_jobs=-1)
        knn_adj_mat_fix = knn_graph.toarray() + knn_graph.toarray().T
        knn_edge_index = np.argwhere(knn_adj_mat_fix > 0)
        edge_index_filename = os.path.join(output_tile_dir, "edge_index.txt")
        np.savetxt(edge_index_filename, knn_edge_index,
                   delimiter='\t', fmt='%i')

    return skipped_points_count


def generate_topology_structures(data_source, input_folder, knn_k, output_folder_root):
    """
    Generate KNN topology structures for each region and partition.

    Args:
        data_source (str): Root directory of the data source.
        input_folder (str): Input directory containing partitioned data.
        knn_k (int): Number of neighbors in KNN.
        output_folder_root (str): Root directory for output.
    """
    region_name_list = load_region_names(data_source)
    output_folder = prepare_output_folder(output_folder_root, knn_k, reset=True)

    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
          "Constructing topology structures of KNN graphs...")

    for region_name in region_name_list["slide"]:
        for partition_type in ['Square', 'Vertical', 'Horizontal']:
            tile_root = os.path.join(input_folder, region_name, partition_type)
            tile_name_list = pd.read_csv(os.path.join(
                tile_root, "tile_name_list.txt"), sep="\t", header=None, names=["tile_name"])
            skipped_points_count = process_each_tile(
                tile_name_list, tile_root, knn_k, output_folder, region_name, partition_type)
            print(
                f"Number of points skipped in the last tile for {region_name} in {partition_type} partition: {skipped_points_count}")
            if skipped_points_count > 0:
                skipped_points_filename = os.path.join(
                    output_folder, region_name, partition_type, "skipped_points.txt")
                with open(skipped_points_filename, 'w') as f:
                    f.write(
                        f"Number of points skipped in the last tile: {skipped_points_count}\n")

    print("All topology structures have been generated!")
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


def generate_node_attribute_matrices(data_source, input_folder, output_folder_root, knn_k):
    """
    Generate node attribute matrices for each KNN graph.

    Args:
        data_source (str): Root directory of the data source.
        input_folder (str): Input directory containing partitioned data.
        output_folder_root (str): Root directory for output.
        knn_k (int): Number of neighbors in KNN.
    """
    region_name_list = load_region_names(data_source)
    output_folder = prepare_output_folder(output_folder_root, knn_k)

    print("Generating node attribute matrices of KNN graphs...")

    for region_name in region_name_list["slide"]:
        for partition_type in ['Square', 'Vertical', 'Horizontal']:
            tile_root = os.path.join(input_folder, region_name, partition_type)
            tile_name_list = pd.read_csv(os.path.join(
                tile_root, "tile_name_list.txt"), sep="\t", header=None, names=["tile_name"])

            cell_type_filename = os.path.join(
                tile_root, "cell_type_encoding.txt")
            n_cell_types = sum(1 for _ in open(cell_type_filename))

            for tile_name in tile_name_list["tile_name"]:
                cell_type_filename = os.path.join(
                    tile_root, f"tile_{tile_name}_cell_type.txt")
                cell_type_label = pd.read_csv(
                    cell_type_filename, sep="\t", header=None, names=["cell_type"])

                node_attr_matrix = np.zeros(
                    (len(cell_type_label), n_cell_types))
                for cell_ind in range(len(cell_type_label)):
                    type_index = cell_type_label["cell_type"][cell_ind]
                    node_attr_matrix[cell_ind, type_index] = 1

                filename1 = os.path.join(
                    output_folder, region_name, partition_type, f"tile_{tile_name}/node_attribute.txt")
                np.savetxt(filename1, node_attr_matrix,
                           delimiter='\t', fmt='%i')

    print("All node attribute matrices have been generated!")
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


class SpatialOmicsImageDataset(InMemoryDataset):
    """
    Custom dataset class for spatial omics data.
    """

    def __init__(self, root, transform=None, pre_transform=None, data_list=None):
        """
        Initialize the dataset, processing if data_list is provided.

        Args:
            root (str): Root directory for the dataset.
            transform (callable, optional): Transformation to apply to the data.
            pre_transform (callable, optional): Transformation to apply before processing.
            data_list (list, optional): List of data objects to initialize the dataset with.
        """
        self.data_list = data_list if data_list is not None else []
        super(SpatialOmicsImageDataset, self).__init__(
            root, transform, pre_transform)
        if self.data_list:
            self.process()
        else:
            self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        """
        List of raw file names. (Not used in this implementation)

        Returns:
            list: Empty list.
        """
        return []

    @property
    def processed_file_names(self):
        """
        List of processed file names.

        Returns:
            list: List containing the processed dataset file name.
        """
        return ['SpatialOmicsImageDataset.pt']

    def download(self):
        """No download functionality needed."""
        pass

    def process(self):
        """
        Process the dataset and save the results.
        """
        data, slices = self.collate(self.data_list)
        torch.save((data, slices), self.processed_paths[0])


def transform_graph_data_structure(data_source, input_folder, output_folder_root, knn_k):
    """
    Transform the graph data structure for each region and partition.

    Args:
        data_source (str): Root directory of the data source.
        input_folder (str): Input directory containing partitioned data.
        output_folder_root (str): Root directory for output.
        knn_k (int): Number of neighbors in KNN.
    """
    region_name_list = load_region_names(data_source)
    output_folder = prepare_output_folder(output_folder_root, knn_k)

    print("Start graph data structure transformation...")

    for region_name in region_name_list["slide"]:
        for partition_type in ['Square', 'Vertical', 'Horizontal']:
            tile_root = os.path.join(input_folder, region_name, partition_type)
            tile_name_list = pd.read_csv(os.path.join(
                tile_root, "tile_name_list.txt"), sep="\t", header=None, names=["tile_name"])

            for tile_name in tile_name_list["tile_name"]:
                output_tile_root = os.path.join(
                    output_folder, region_name, partition_type, f"tile_{tile_name}")

                edge_index_filename = os.path.join(
                    output_tile_root, "edge_index.txt")
                if os.path.exists(edge_index_filename):
                    edge_index = torch.from_numpy(np.loadtxt(
                        edge_index_filename, dtype='int64', delimiter="\t").T)

                    node_attribute_filename = os.path.join(
                        output_tile_root, "node_attribute.txt")
                    x = torch.from_numpy(np.loadtxt(
                        node_attribute_filename, dtype='float32', delimiter="\t"))

                    data = Data(x=x, edge_index=edge_index.contiguous())
                    SpatialOmicsImageDataset(
                        output_tile_root, transform=T.ToDense(), data_list=[data])

    print("Graph data structure transformation complete!")
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == "__main__":
    args = parse_arguments()

    generate_topology_structures(
        args.data_source, args.input_folder, args.knn_k, args.output_folder)
    generate_node_attribute_matrices(
        args.data_source, args.input_folder, args.output_folder, args.knn_k)
    transform_graph_data_structure(
        args.data_source, args.input_folder, args.output_folder, args.knn_k)
