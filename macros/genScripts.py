# Create scripts for input into Geant4 with CHPC jop Array
import subprocess
import os

def main () :
  energies = ["1.00e-7", "4.14e-7", "8.76e-7", "1.86e-6", "5.04e-6", "1.07e-5", "3.73e-5", "1.01e-4", "2.14e-4", "4.54e-4", "1.58e-3", "3.35e-3", "7.1e-3", "1.5e-2", "2.19e-2", "2.42e-2", "2.61e-2", "3.18e-2", "4.09e-2", "6.74e-2", "1.11e-1", "1.83e-1", "2.97e-1", "3.69e-1", "4.98e-1", "6.08e-1", "7.43e-1", "8.23e-1", "1.0", "1.35", "1.65", "1.92", "2.23", "2.35", "2.37", "2.47", "2.73", "3.01", "3.68", "4.97", "6.07", "7.41", "8.61", "10.0", "12.2", "14.2", "17.3"]

  # Get the starter file that has the default input for the script:
  readFile = open('bf3.mac', 'r')
  for i in range(len(energies)) :
    readFile = open('bf3.mac', 'r')
    # Open the file with the name I want to write to:
    indexName = str(i)
    writeFile = open('group-'+indexName+'.mac', 'w')
    writeFile.write("/RunAction/FileName group-" + indexName +"\n")
    for line in readFile.readlines() :
      if line.split(' ')[0] == "/gps/energy" :
        writeFile.write("/gps/energy " + energies[i] + " MeV\n")
      else :
        writeFile.write(line)
    readFile.close()  
    writeFile.close()

if __name__ == '__main__' :
  main()
