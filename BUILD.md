# Instructions for building binaries
i.e. `.exe` files (if you're using Windows)

I don't have a Mac so I can't make a MacOS executable. If you need one, reach out
and we can arrange something

## Linux

1. **Create a 'hook' file:**
```python
# /hooks/hook_.py
from PyInstaller.utils.hooks import copy_metadata
datas = []
datas += copy_metadata('pyopenms')
```

2. **Run pyinstaller:**
```bash
pyinstaller --onefile --additional-hooks-dir ./hooks waters_mzml_tools.py --clean
```

3. **Modify `datas` in the `.spec` file:**
```python
# ...
# datas=[],
#           Change to..
datas=[
    (".venv/lib/python3.10/site-packages/pyopenms", "./pyopenms/"),
],
```

4. **Run pyinstaller again:**
```bash
pyinstaller waters_mzml_tools.spec --clean
```

---

## Windows
1. Set up virtual environment using `venv` and `pip install -r requirements.txt`
2. Go to `<Project_dir>\venv-3.10\Lib\site-packages\PyQt5\Qt5\bin` and find:
    - `Qt5Core.dll`
    - `Qt5Network.dll`
3. Copy those files to `<Project_dir>\venv-3.10\Lib\site-packages\pyopenms`
    - Feel free to rename the old `.dll`'s instead of over-writing them
4. Now follow the same steps as in the *Linux* section 
