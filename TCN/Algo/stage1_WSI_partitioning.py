import os
import math
import shutil
import datetime
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict


def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Partition and save sample data based on different partition types.")
    parser.add_argument("--data_root", type=str, default="./Processed_data/",
                        help="Root directory containing sample data files.")
    parser.add_argument("--min_size", type=int, default=100,
                        help="Minimum number of cells allowed in a tile.")
    parser.add_argument("--max_size", type=int, default=8000,
                        help="Maximum number of cells allowed in a tile.")
    parser.add_argument("--output_folder", type=str, default="./stage0_output/",
                        help="Output folder to save the partitioned data.")
    return parser.parse_args()


def load_data(file_path, column_names):
    """
    Load data from a file into a DataFrame.

    Args:
        file_path (str): Path to the file.
        column_names (list): List of column names for the DataFrame.

    Returns:
        pd.DataFrame: Loaded data.
    """
    return pd.read_csv(file_path, delimiter="\t", header=None, names=column_names)


def calculate_tile_dimensions(data_shape, max_size, x_range, y_range, partition_type):
    """
    Calculate the dimensions of tiles based on the partition type.

    Args:
        data_shape (tuple): Shape of the data.
        max_size (int): Maximum number of cells allowed in a tile.
        x_range (tuple): Range of x-coordinates.
        y_range (tuple): Range of y-coordinates.
        partition_type (str): Type of partition (Square, Vertical, Horizontal).

    Returns:
        tuple: Number of tiles in x and y directions, and the size of each tile.
    """
    if partition_type == 'Square':
        n_tiles_x = n_tiles_y = int(math.sqrt(data_shape[0] / max_size)) + 1
    elif partition_type == 'Vertical':
        n_tiles_x, n_tiles_y = int(data_shape[0] / max_size) + 1, 1
    elif partition_type == 'Horizontal':
        n_tiles_x, n_tiles_y = 1, int(data_shape[0] / max_size) + 1
    cell_size_x = (x_range[1] - x_range[0]) / n_tiles_x
    cell_size_y = (y_range[1] - y_range[0]) / n_tiles_y
    return n_tiles_x, n_tiles_y, cell_size_x, cell_size_y


def split_tile_recursively(tile_coords, tile_info, max_size, base_cell_str=""):
    """
    Recursively split tiles that exceed the maximum size.

    Args:
        tile_coords (np.ndarray): Coordinates of cells in the tile.
        tile_info (list): Cell type information for the tile.
        max_size (int): Maximum number of cells allowed in a tile.
        base_cell_str (str): Base string for naming split tiles.

    Returns:
        list: List of split tiles with coordinates, info, and names.
    """
    split_tiles = []
    n_splits = math.ceil(len(tile_coords) / max_size)
    for i in range(n_splits):
        start_idx, end_idx = i * max_size, (i + 1) * max_size
        sub_tile_coords, sub_tile_info = tile_coords[start_idx:
                                                     end_idx], tile_info[start_idx:end_idx]
        sub_tile_str = f"{base_cell_str}_{i}"
        if len(sub_tile_coords) > max_size:
            sub_split_tiles = split_tile_recursively(
                sub_tile_coords, sub_tile_info, max_size, sub_tile_str)
            split_tiles.extend(sub_split_tiles)
        else:
            split_tiles.append((sub_tile_coords, sub_tile_info, sub_tile_str))
    return split_tiles


def merge_small_tiles(tiles, info_tiles, min_size):
    """
    Merge tiles that are smaller than the minimum size.

    Args:
        tiles (dict): Dictionary of tile coordinates.
        info_tiles (dict): Dictionary of tile cell type information.
        min_size (int): Minimum number of cells allowed in a tile.
    """
    small_tiles = [cell for cell,
                   coords in tiles.items() if len(coords) < min_size]

    for cell in small_tiles:
        index = list(tiles.keys()).index(cell)
        merged = False
        neighbors = get_next_neighbors(
            index, tiles) if index == 0 else get_previous_neighbors(index)
        for neighbor_idx in neighbors:
            neighbor = list(tiles.keys())[neighbor_idx]
            if neighbor in tiles:
                tiles[neighbor].extend(tiles[cell])
                info_tiles[neighbor].extend(info_tiles[cell])
                merged = True
                break
        if merged:
            del tiles[cell]
            del info_tiles[cell]


def get_previous_neighbors(index):
    """
    Get indices of previous neighboring tiles.

    Args:
        index (int): Index of the current tile.

    Returns:
        list: List of indices of previous neighbors.
    """
    return [index - i for i in range(1, index + 1)]


def get_next_neighbors(index, tiles):
    """
    Get indices of next neighboring tiles.

    Args:
        index (int): Index of the current tile.
        tiles (dict): Dictionary of tile coordinates.

    Returns:
        list: List of indices of next neighbors.
    """
    return [index + i for i in range(1, len(tiles) - index)]


def save_partitioned_tiles(output_dir, all_tiles, all_info_tiles):
    """
    Save partitioned tiles and their information to files.

    Args:
        output_dir (str): Directory to save the tiles.
        all_tiles (dict): Dictionary of all tile coordinates.
        all_info_tiles (dict): Dictionary of all tile cell type information.

    Returns:
        tuple: List of tile names, tile order, and tile size information.
    """
    os.makedirs(output_dir, exist_ok=True)
    tile_names, tile_order, tile_size_info = [], [], []

    for cell_str, coords in all_tiles.items():
        coords_sorted = np.array(coords)
        coords_sorted = coords_sorted[np.lexsort(
            (coords_sorted[:, 1], coords_sorted[:, 0]))]

        info_sorted = [info for _, info in sorted(
            zip(coords, all_info_tiles[cell_str]), key=lambda pair: (pair[0][0], pair[0][1]))]

        coords_file = os.path.join(
            output_dir, f"tile_{cell_str}_coordinate.txt")
        info_file = os.path.join(output_dir, f"tile_{cell_str}_cell_type.txt")

        np.savetxt(coords_file, coords_sorted, fmt='%f', delimiter="\t")
        np.savetxt(info_file, info_sorted, fmt='%d')

        tile_names.append(cell_str)
        tile_order.append(cell_str)
        tile_size_info.append(f"{cell_str}\t{len(coords)}")

    return tile_names, tile_order, tile_size_info


def compare_dataframes(original_coords_sorted, resequenced_coords_sorted):
    """
    Compare original and resequenced coordinates to ensure consistency.

    Args:
        original_coords_sorted (pd.DataFrame): Sorted original coordinates.
        resequenced_coords_sorted (pd.DataFrame): Sorted resequenced coordinates.
    """
    def dataframes_are_close(df1, df2, rtol=1e-05, atol=10000):
        return np.allclose(df1.to_numpy(), df2.to_numpy(), rtol=rtol, atol=atol)

    coords_match = dataframes_are_close(
        original_coords_sorted[['x', 'y']], resequenced_coords_sorted[['x', 'y']])
    if not coords_match:
        differences = original_coords_sorted.compare(
            resequenced_coords_sorted, keep_shape=True, keep_equal=False)
        print("Differences in coordinates detected:")
        print(differences)
    assert coords_match, "Mismatch in coordinates content detected."


def partition_and_save(data_root, region_name_list, partition_type, max_size, min_size, output_folder):
    """
    Main function to partition and save data into tiles based on the partition type.

    Args:
        data_root (str): Root directory containing data files.
        region_name_list (pd.DataFrame): List of sample names.
        partition_type (str): Type of partition (Square, Vertical, Horizontal).
        max_size (int): Maximum number of cells allowed in a tile.
        min_size (int): Minimum number of cells allowed in a tile.
        output_folder (str): Directory to save the partitioned data.
    """
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    print(f"Partitioning and saving input for {partition_type} partition...")

    for graph_index, region_name in enumerate(region_name_list.slide):
        print(f"Processing sample-{graph_index}")

        coordinates_file_path = os.path.join(
            data_root, f"{region_name}_coordinate.txt")
        info_file_path = os.path.join(
            data_root, f"{region_name}_cell_type.txt")

        data = load_data(coordinates_file_path, ['x', 'y'])
        info = load_data(info_file_path, ['info'])
        original_data = data.copy()

        x_min, x_max = data['x'].min(), data['x'].max()
        y_min, y_max = data['y'].min(), data['y'].max()

        n_tiles_x, n_tiles_y, cell_size_x, cell_size_y = calculate_tile_dimensions(
            data.shape, max_size, (x_min, x_max), (y_min, y_max), partition_type)

        coordinates = data[['x', 'y']].to_numpy()
        info_array = info['info'].to_numpy()
        info_map = {text: idx for idx,
                    text in enumerate(np.unique(info_array))}
        info_int_array = np.vectorize(info_map.get)(info_array)

        x_indices = np.floor(
            (coordinates[:, 0] - x_min) / cell_size_x).astype(int)
        y_indices = np.floor(
            (coordinates[:, 1] - y_min) / cell_size_y).astype(int)

        cell_ids = np.stack((x_indices, y_indices), axis=-1)
        tiles, info_tiles = defaultdict(list), defaultdict(list)

        # Assign cells to tiles
        for idx, cell in enumerate(cell_ids):
            cell_tuple = tuple(cell)
            tiles[cell_tuple].append(coordinates[idx])
            info_tiles[cell_tuple].append(info_int_array[idx])

        all_tiles, all_info_tiles = defaultdict(list), defaultdict(list)
        # Split large tiles
        for cell, coords in tiles.items():
            if len(coords) > max_size:
                split_tiles = split_tile_recursively(
                    coords, info_tiles[cell], max_size, f"{cell[0]}_{cell[1]}")
                for split_coords, split_info, cell_str in split_tiles:
                    all_tiles[cell_str] = split_coords
                    all_info_tiles[cell_str] = split_info
            else:
                cell_str = f"{cell[0]}_{cell[1]}"
                all_tiles[cell_str] = coords
                all_info_tiles[cell_str] = info_tiles[cell]

        # Merge small tiles
        merge_small_tiles(all_tiles, all_info_tiles, min_size)

        output_dir = os.path.join(output_folder, region_name, partition_type)
        tile_names, tile_order, tile_size_info = save_partitioned_tiles(
            output_dir, all_tiles, all_info_tiles)

        tile_name_list_file = os.path.join(output_dir, 'tile_name_list.txt')
        with open(tile_name_list_file, 'w') as f:
            f.write("\n".join(tile_names))

        tile_size_info_file = os.path.join(output_dir, 'tile_size_info.txt')
        with open(tile_size_info_file, 'w') as f:
            f.write("\n".join(tile_size_info))

        info_map_file = os.path.join(output_dir, 'cell_type_encoding.txt')
        with open(info_map_file, 'w') as f:
            for text, idx in info_map.items():
                f.write(f"{idx} {text}\n")

        print(f"Partitioning and saving completed for sample: {graph_index}")

        # Resequence tiles into the original sequence
        resequenced_data = []
        for cell_str in tile_order:
            coords_file = os.path.join(
                output_dir, f"tile_{cell_str}_coordinate.txt")
            info_file = os.path.join(
                output_dir, f"tile_{cell_str}_cell_type.txt")

            if os.path.exists(coords_file) and os.path.exists(info_file):
                coords = np.loadtxt(coords_file, delimiter="\t")
                info = np.loadtxt(info_file, delimiter="\t")

                coords_df = pd.DataFrame(coords, columns=['x', 'y'])
                info_df = pd.DataFrame(info, columns=['info'])
                combined_df = pd.concat([coords_df, info_df], axis=1)
                resequenced_data.append(combined_df)

        final_resequenced_data = pd.concat(resequenced_data, ignore_index=True)
        original_coords_sorted = original_data[['x', 'y']].sort_values(
            by=['x', 'y']).reset_index(drop=True)
        resequenced_coords_sorted = final_resequenced_data[[
            'x', 'y']].sort_values(by=['x', 'y']).reset_index(drop=True)

        compare_dataframes(original_coords_sorted, resequenced_coords_sorted)
        print(
            f"Total coordinates match for sample {graph_index} after resequencing.")

    print(
        f"All samples processed successfully for {partition_type} partition.")


if __name__ == "__main__":
    args = parse_arguments()

    region_filename = os.path.join(args.data_root, "sample_name_list.txt")
    region_name_list = pd.read_csv(
        region_filename, sep="\t", header=None, names=["slide"])

    if os.path.exists(args.output_folder):
        shutil.rmtree(args.output_folder)
    os.makedirs(args.output_folder)

    partition_types = ['Square', 'Vertical', 'Horizontal']
    for partition_type in partition_types:
        partition_and_save(args.data_root, region_name_list, partition_type,
                           args.max_size, args.min_size, args.output_folder)
