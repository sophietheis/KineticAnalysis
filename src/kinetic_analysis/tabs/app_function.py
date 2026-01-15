import io
import os
import base64

import tkinter as tk
from tkinter import filedialog

from ..utils.utils import read_csv_file


def upload_csv(contents, app, name="csv_to_analyse"):
    if contents is None:
        return None, ""

    # Decode and read the CSV content
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        df = read_csv_file(io.StringIO(decoded.decode('utf-8')))

    except Exception as e:
        return None, f"Failed to parse CSV: {str(e)}"

    # Save to app data
    app.data[name] = df

    return df, f"Success to parse CSV"


def browse_directory(n_clicks, col_name, app):
    if n_clicks:
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        if folder_selected:
            print(folder_selected)
        else:
            print(None)
        app.data[col_name] = folder_selected
        return f"Directory chosen: {app.data[col_name]}"


def list_csv_files(directory, col_name, app):
    """
    List csv files inside the directory.
    List of files is stored in col_name.
    :param directory:
    :type directory:
    :param col_name:
    :type col_name:
    :param app:
    :type app:
    :return:
    :rtype:
    """
    if app.data[col_name]:
        app.data['csv_files'] = [
            {'label': file, 'value': file}
            for file in os.listdir(app.data[col_name]) if
            file.endswith('.csv')
        ]
        return app.data['csv_files']
    return []
