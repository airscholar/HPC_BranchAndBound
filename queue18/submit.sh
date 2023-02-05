search_dir=.
for entry in "$search_dir"/*
do
  #if entry name ends with .sub
  if [[ $entry == *.sub ]]; then
    #remove ./ from entry
    entry=${entry:2}
    echo "Submitting $entry"
    qsub $entry
  fi

#  echo "qsub $entry"
done