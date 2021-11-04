mkdir all_gifs
for d in ./*/ ; do 
	cd "$d"
	VALUE=$(basename $d)
	rm files.txt
	rm "$VALUE".gif
	rm concat.avi 
	for f in *.avi; do 
		echo "file $f" >> files.txt
	done
        ffmpeg -f concat -safe 0 -i files.txt -c copy concat.avi	
	ffmpeg -i concat.avi -vf "fps=30,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 "$VALUE".gif
	cp "$VALUE".gif ../all_gifs/
	cd ..
done

