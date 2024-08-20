# run this command in the plugin directory
plugins=$(cat ../build/plugins.txt)
echo $plugins
~/software/checkheaders-1.0.1/checkheaders $plugins  > tmp_checkheader_progress 2> tmp_checkheader_errors
grep -v "Header not found" tmp_checkheader_errors | grep -v "test" > tmp_checkheader_errors2
