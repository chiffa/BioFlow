def _is_int(_obj):
    """
    Checks if an object is an int with a try-except loop

    :param _obj:
    :return:
    """
    try:
        int(_obj)
    except TypeError or ValueError as e:
        return False
    else:
        return True