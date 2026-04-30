rm -rf $1
mkdir $1
cp B-Spline-IPC/config.yaml $1
mkdir $1/models
mkdir $1/quad_data
mkdir $1/seam
cp -r B-Spline-IPC/models/* $1/models
cp -r B-Spline-IPC/quad_data/* $1/quad_data
cp -r B-Spline-IPC/seam/* $1/seam
cp bin/Release-linux-x86_64/B-Spline-IPC/* $1
cd $1
./B-Spline-IPC
pause