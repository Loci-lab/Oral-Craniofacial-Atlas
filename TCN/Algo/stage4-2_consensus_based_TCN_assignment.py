import os
import pandas as pd
import numpy as np
import argparse
import shutil


def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process cell type encoding and generate result table.")

    parser.add_argument("--image_name", type=str, default="combined_all_PCF",
                        help="Name of the image data to be processed.")
    parser.add_argument("--stage1_root", type=str, default="./Results/Intermediate_ouputs/stage1_output/",
                        help="Root directory for stage 1 outputs.")
    parser.add_argument("--last_stage_root", type=str, default="./Results/Intermediate_ouputs/stage4-1_output/",
                        help="Root directory for last stage outputs.")
    parser.add_argument("--this_stage_root", type=str, default="./Results/Intermediate_ouputs/stage4-2_output/",
                        help="Root directory for this stage outputs.")
    parser.add_argument("--load_config", type=str, default="k-69_TCN-20_lr-0.003/",
                        help="Configuration of the model to load.")
    parser.add_argument("--epoch_to_use", type=int, default=480,
                        help="Epoch to use for the model.")

    return parser.parse_args()


def read_and_merge_files(filenames, column_names):
    """
    Read and merge multiple CSV files into a single DataFrame.

    Args:
        filenames (list): List of file paths to read.
        column_names (list): Column names for the resulting DataFrame.

    Returns:
        pd.DataFrame: Merged DataFrame containing data from all files.
    """
    dataframes = [pd.read_csv(
        file, sep="\t", header=None, names=column_names) for file in filenames]
    return pd.concat(dataframes, ignore_index=True)


def read_data(partition_type, folder):
    """
    Read and sort data from a CSV file based on partition type.

    Args:
        partition_type (str): The partition type (Square, Vertical, Horizontal).
        folder (str): The folder containing the result files.

    Returns:
        pd.DataFrame: Sorted DataFrame.
    """
    filename = os.path.join(folder, f"{partition_type}_result_table.csv")
    df = pd.read_csv(filename, sep=",", header=0, names=["x_coordinate", "y_coordinate", "Cell_Type", "TCN_by_{}".format(partition_type)], dtype={
                     "x_coordinate": np.float32, "y_coordinate": np.float32, "Cell_Type": np.int32, f"TCN_by_{partition_type}": np.int32})
    return df.sort_values(by=['x_coordinate', 'y_coordinate']).reset_index(drop=True)


def majority_vote(row):
    """
    Determine the majority vote among the TCN labels from different partitions.

    Args:
        row (pd.Series): Row containing TCN labels from Square, Vertical, and Horizontal partitions.

    Returns:
        int: The majority vote label.
    """
    labels = [row['TCN_by_Square'],
              row['TCN_by_Horizontal'], row['TCN_by_Vertical']]
    return int(max(set(labels), key=labels.count))


def main(args):
    # Define output folders
    last_stage_output_folder = os.path.join(
        args.last_stage_root, args.load_config, args.image_name, f"checkpoint_epoch_{args.epoch_to_use}")
    this_stage_output_folder = os.path.join(
        args.this_stage_root, args.load_config, args.image_name, f"checkpoint_epoch_{args.epoch_to_use}")

    if os.path.exists(this_stage_output_folder):
        shutil.rmtree(this_stage_output_folder)
    os.makedirs(this_stage_output_folder, exist_ok=True)

    for partition_type in ["Square", "Vertical", "Horizontal"]:
        # Prepare file paths for coordinates and cell types
        tile_root = os.path.join(
            args.stage1_root, args.image_name, partition_type)
        tile_name_list = pd.read_csv(os.path.join(
            tile_root, "tile_name_list.txt"), sep="\t", header=None, names=["tile_name"])

        coordinate_filenames = [os.path.join(
            tile_root, f"tile_{tile_name}_coordinate.txt") for tile_name in tile_name_list["tile_name"]]
        cell_type_filenames = [os.path.join(
            tile_root, f"tile_{tile_name}_cell_type.txt") for tile_name in tile_name_list["tile_name"]]

        x_y_coordinates = read_and_merge_files(
            coordinate_filenames, ["x_coordinate", "y_coordinate"])
        cell_types = read_and_merge_files(cell_type_filenames, ["cell_type"])

        # Load node mask and cluster assignment matrix
        node_mask = pd.read_csv(os.path.join(
            last_stage_output_folder, f"{partition_type}_node_mask.csv"), header=None)
        nonzero_ind = node_mask[node_mask.iloc[:, 0] == 1].index

        soft_clust_file = os.path.join(
            last_stage_output_folder, f"{partition_type}_TCN_assign_matrix.csv")
        clust_matrix = pd.read_csv(
            soft_clust_file, header=None).iloc[nonzero_ind, :]
        hard_clust_label = clust_matrix.idxmax(axis=1).values
        hard_clust_label = pd.DataFrame(hard_clust_label, columns=["TCN"])

        # Create result table and save
        result_table = pd.concat(
            [x_y_coordinates, cell_types, hard_clust_label], axis=1)
        result_table.to_csv(os.path.join(this_stage_output_folder,
                            f"{partition_type}_result_table.csv"), index=False, header=True)

    # Read data from all partitions
    square_df = read_data("Square", this_stage_output_folder)
    horizontal_df = read_data("Horizontal", this_stage_output_folder)
    vertical_df = read_data("Vertical", this_stage_output_folder)

    # Ensure data consistency across partitions
    assert square_df.shape == horizontal_df.shape == vertical_df.shape

    # Merge data from all partitions and apply majority vote
    df = pd.merge(square_df, horizontal_df, on=[
                  "x_coordinate", "y_coordinate", "Cell_Type"], how="inner")
    df = pd.merge(df, vertical_df, on=[
                  "x_coordinate", "y_coordinate", "Cell_Type"], how="inner")
    
    condition = (df['TCN_by_Square'] != df['TCN_by_Horizontal']) & \
                (df['TCN_by_Square'] != df['TCN_by_Vertical']) & \
                (df['TCN_by_Horizontal'] != df['TCN_by_Vertical'])

    # Removing those rows
    df = df[~condition]

    df['TCN'] = df.apply(majority_vote, axis=1)

    # Load and sort unique cell types
    unique_cell_type_df = pd.read_csv(os.path.join(args.stage1_root, args.image_name, "Square",
                                      "cell_type_encoding.txt"), sep="\t", header=None, names=["cell_type_encoding"])
    cell_type_encoding = unique_cell_type_df['cell_type_encoding'].values.tolist(
    )

    sorted_unique_labels = [f"{int(item.split()[0])} {item.split(' ', 1)[1]}" for item in sorted(
        cell_type_encoding, key=lambda x: int(x.split()[0]))]

    # Map cell types to their string labels
    df_output = df.loc[:, ["x_coordinate", "y_coordinate", "Cell_Type", "TCN"]]
    df_output["Cell_Type"] = df_output["Cell_Type"].apply(
        lambda x: sorted_unique_labels[int(x)])

    # Save final result table
    df_output.to_csv(os.path.join(this_stage_output_folder,
                     "final_result_table.csv"), index=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
