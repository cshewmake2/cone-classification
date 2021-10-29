mkdir all_gifs
for d in ./*/ ; do 
	cd "$d"
	rm files.txt
	rm concat.avi 
	for f in *.avi; do 
		echo "file $f" >> files.txt
	done
	VALUE=$(basename $d)
        ffmpeg -f concat -safe 0 -i files.txt -c copy concat.avi	
	ffmpeg  -t 5 -i concat.avi -vf "fps=30,scale=320:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 "$VALUE".gif
	cp "$VALUE".gif ../all_gifs/
	cd ..
done

