import scipy.io

def read_mat_file(filename):
    """ Load the .mat file into a dictionary
    Example:
    filename = 'example.mat'
    variable_names, variable_types, variable_values = read_mat_file(filename)
    print(variable_names)
    print(variable_types)
    print(variable_values)
    """
    mat_data = scipy.io.loadmat(filename)

    # Get a list of variable names
    variable_names = [var for var in mat_data if not var.startswith('__')]

    # Get a list of variable types and values
    variable_types = [type(mat_data[var]) for var in variable_names]
    variable_values = [mat_data[var] for var in variable_names]

    return variable_names, variable_types, variable_values