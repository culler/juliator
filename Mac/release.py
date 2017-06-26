import os, sys, re
from subprocess import Popen, call, PIPE
from math import ceil

def get_tk_ver(python):
    """
    Figure out which version of Tk is used by this python.
    """
    out, errors = Popen([python, "-c", "import _tkinter; print(_tkinter.TK_VERSION)"],
                    stdout = PIPE).communicate()
    return out.strip()

def freshen_app(python):
    """
    Pull and build the module
    """
    os.chdir("../")
    call(["hg", "pull"])
    call(["hg", "up"])
    call([python, "setup.py", "install"])
    os.chdir("Mac")

def build_app(python):
    """
    Build the standalone app bundle.
    """
    call([python, "setup.py", "py2app"])

def cleanup_app(python):
    """
    Tidy things up.
    """
    # Add a symlink so that py2app 0.13 can find the Tcl Scripts directory
    tk_ver = get_tk_ver(python)
    libdir = "dist/Juliator.app/Contents/lib/"
    scriptdir = "../Frameworks/Tk.Framework/Versions/%s"%tk_ver
    os.mkdir(libdir)
    os.symlink(scriptdir, libdir + "tk%s"%tk_ver)

def package_app(dmg_name):
    """
    Create a disk image containing the app, with a nice background and
    a symlink to the Applications folder.
    """
    image_dir = "disk_images"
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    mount_name = os.path.join("/Volumes", dmg_name)
    dmg_path = os.path.join(image_dir, dmg_name + ".dmg")
    temp_path = os.path.join(image_dir, dmg_name + "-tmp.dmg")
    # Make sure the dmg isn't currently mounted, or this won't work.  
    while os.path.exists(mount_name):
        print("Trying to eject " + mount_name)
        os.system("hdiutil detach " + mount_name)
    # Remove old dmgs if they exist.
    if os.path.exists(dmg_path):
        os.remove(dmg_path)
    if os.path.exists(temp_path):
        os.remove(temp_path)
    # Add symlink to /Applications if not there.
    if not os.path.exists("dist/Applications"):
        os.symlink("/Applications/", "dist/Applications")

    # Copy over the background and .DS_Store file.
    call(["rm", "-rf", "dist/.background.png"])
    call(["cp", "background.png", "dist/.background.png"])
    call(["cp", "dotDS_Store", "dist/.DS_Store"])
        
    # Figure out the needed size.
    raw_size, errors = Popen(["du", "-sh", "dist"], stdout=PIPE).communicate()
    size, units = re.search("([0-9.]+)([KMG])", raw_size).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units
    # Run hdiutil to build the dmg file.:
    call(["hdiutil", "makehybrid", "-hfs", "-hfs-volume-name", dmg_name,
        "-hfs-openfolder", "dist", "dist", "-o", temp_path])
    call(["hdiutil", "convert", "-format", "UDZO", temp_path, "-o", dmg_path])
    os.remove(temp_path)
    # Delete the symlink to /Applications or egg_info will be glacial on newer setuptools.
    os.remove("dist/Applications")

def do_release(python, dmg_name):
    freshen_app(python)
    build_app(python)
    cleanup_app(python)
    package_app(dmg_name)

do_release('python3', 'Juliator')

