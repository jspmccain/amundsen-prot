### Folder to look through the mzML files
folder='../data/mzML-converted'

### Loop through this folder and run dinosaur on every mzML file
for file in "$folder"/*.mzML; do
    echo $file
    java -jar ../../Dinosaur-1.2.0.free.jar --verbose --profiling --concurrency=4 --outDir=../data/dino-converted/ $file
done
