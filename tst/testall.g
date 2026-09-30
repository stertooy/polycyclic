LoadPackage("polycyclic");
Print("DEBUG VARIABLES ON: ",CHECK_IGS@Polycyclic,"\n");
TestDirectory(DirectoriesPackageLibrary("polycyclic", "tst"), rec(exitGAP := true));
FORCE_QUIT_GAP(1);
