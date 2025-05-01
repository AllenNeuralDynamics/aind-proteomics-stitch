dir = getDirectory("plugins");
print("Plugins directory: " + dir);
list = getFileList(dir);
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".class") || endsWith(list[i], ".jar") || endsWith(list[i], ".py") || endsWith(list[i], ".ijm")) {
        print(list[i]);
    }
}