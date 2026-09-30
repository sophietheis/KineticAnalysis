import io
import os
import base64
import logging

import tkinter as tk
from tkinter import filedialog

from ..utils.utils import read_csv_file
from ..session_store import get_session_data, set_session_data

logger = logging.getLogger(__name__)


def upload_csv(contents, session_id, name="csv_to_analyse"):
    if contents is None:
        return None, ""

    # Decode and read the CSV content
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        df = read_csv_file(io.StringIO(decoded.decode('utf-8')))

    except Exception as e:
        logger.exception("Failed to parse uploaded CSV for '%s'", name)
        return None, f"Failed to parse CSV: {str(e)}"

    # Save to this browser session's data only
    set_session_data(session_id, name, df)

    return df, f"Success to parse CSV"


def browse_directory(n_clicks, col_name, session_id):
    if n_clicks:
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        folder_selected = filedialog.askdirectory()
        root.destroy()
        logger.debug("Directory selected for '%s': %s", col_name, folder_selected or None)
        
        set_session_data(session_id, col_name, folder_selected)
        return f"Directory chosen: {folder_selected}"


def list_csv_files(directory, col_name, session_id):
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
    directory_value = get_session_data(session_id, col_name)
    if directory_value:
        csv_files = [
            {'label': file, 'value': file}
            for file in os.listdir(directory_value) if
            file.endswith('.csv')
        ]
        return csv_files
    return []
