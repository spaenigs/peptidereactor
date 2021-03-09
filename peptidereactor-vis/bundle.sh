##!/usr/bin/env sh

ROOT_DIR=data/temp/peptidereactor-public/

echo "Copying..."

# create dirs in root dir
for dir_path in $(find data/*/vis -type d)
do
  mkdir -p $ROOT_DIR$dir_path
done

# copy files to target dir
for file_path in $(find data/*/vis -type f)
do
  cp $file_path $ROOT_DIR$(dirname $file_path)
done

# copy misc dirs
cp -r peptidereactor-vis/ $ROOT_DIR
cp -r nodes/ $ROOT_DIR

# create README file
README_FILE=$ROOT_DIR/README.md
touch $$README_FILE
cat > $README_FILE <<- EOF
1) \`screen -S peptidereactor # if not running\`

   or:

   \`screen -r peptidereactor\`

2) \`cd peptidereactor-public/\`

3) \`./peptidereactor-vis/run_server\`

4) Ctrl-a d (leave session without quitting)
EOF

echo "Updating..."

SCRIPT_PATH=$ROOT_DIR"peptidereactor-vis/resources/streamlit.py"

FAVICON_URL_LOCALHOST=http://127.0.0.1:8383/favicon.png
FAVICON_URL_PUBLIC=https://peptidereactor.mathematik.uni-marburg.de/favicon.png
sed -i "s%$FAVICON_URL_LOCALHOST%$FAVICON_URL_PUBLIC%g" $ROOT_DIR"peptidereactor-vis/resources/streamlit.py"

LOGO_URL_LOCALHOST=http://127.0.0.1:8383/logo.png
LOGO_URL_PUBLIC=https://peptidereactor.mathematik.uni-marburg.de/logo.png
sed -i "s%$LOGO_URL_LOCALHOST%$LOGO_URL_PUBLIC%g" $ROOT_DIR"peptidereactor-vis/resources/streamlit.py"

echo "Compressing..."

# tar files for upload
tar -czf data/temp/peptidereactor-public.tar.gz -C data/temp/ peptidereactor-public
