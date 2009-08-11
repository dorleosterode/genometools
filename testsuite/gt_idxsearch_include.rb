allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "RandomN.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

maxpairstestfiles=["Duplicate.fna",
                   "Wildcards.fna",
                   "hs5hcmvcg.fna",
                   "humhbb.fna",
                   "mpomtcg.fna",
                   "ychrIII.fna"]

def makegreedyfwdmatcall(queryfile,indexarg,ms)
  prog=""
  if ms
    prog="#{$bin}gt matstat -verify"
  else
    prog="#{$bin}gt uniquesub"
  end
  constantargs="-min 1 -max 20 -query #{queryfile} #{indexarg}"
  return "#{prog} -output querypos #{constantargs}"
end

def checkgreedyfwdmat(queryfile,ms)
  run_test(makegreedyfwdmatcall(queryfile,"-fmi fmi",ms), :maxtime => 1200)
  run "mv #{$last_stdout} tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-esa sfx",ms), :maxtime => 1200)
  run "mv #{$last_stdout} tmp.esa"
  run "diff tmp.esa tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-pck pck",ms), :maxtime => 1200)
  run "mv #{$last_stdout} tmp.pck"
  run "diff tmp.pck tmp.fmi"
end

def checktagerator(queryfile,ms)
  run "#{$bin}gt shredder -minlength 12 -maxlength 15 #{queryfile} | " +
      "#{$bin}gt seqfilter -minlength 12 - | " +
      "sed -e \'s/^>.*/>/\' > patternfile"
  if File.size("patternfile") > 0
    run_test("#{$bin}gt tagerator -rw -cmp -e 0 -esa sfx -q patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -e 1 -esa sfx -q patternfile " + 
             "-withwildcards",:maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -e 2 -esa sfx -q patternfile " +
             "-withwildcards",:maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -esa sfx -q patternfile " +
             " -maxocc 10",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -e 0 -pck pck -q patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -e 1 -pck pck -q patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -e 2 -pck pck -q patternfile",
             :maxtime => 200)
    run_test("#{$bin}gt tagerator -rw -cmp -pck pck -q patternfile " +
             "-maxocc 10",
             :maxtime => 100)
  end
end

def createandcheckgreedyfwdmat(reffile,queryfile)
  run("#{$scriptsdir}/runmkfm.sh #{$bin}/gt 0 . fmi #{reffile}",
      :maxtime => 100)
  run "#{$bin}gt suffixerator -indexname sfx -tis -suf -ssp -dna -v " +
           "-db #{reffile}"
  run("#{$bin}gt packedindex mkindex -tis -ssp -indexname pck -db #{reffile} " +
      "-sprank -dna -pl -bsize 10 -locfreq 32 -dir rev", :maxtime => 100)
  run "#{$bin}gt prebwt -maxdepth 4 -pck pck"
  checkgreedyfwdmat(queryfile,false)
  checkgreedyfwdmat(queryfile,true)
end

def addfilepath(filename)
  if filename == 'Duplicate.fna'
    return "#{$testdata}/#{filename}"
  else
    return "#{$gttestdata}/DNA-mix/Grumbach.fna/#{filename}"
  end
end

def checkmaxpairs(reffile)
  reffilepath=addfilepath(reffile)
  run_test "#{$bin}gt suffixerator -algbds 3 40 120 -db " +
           "#{reffilepath} -indexname sfxidx -dna -suf -tis -lcp -pl"
  resultfile="#{$gttestdata}maxpairs-result/#{reffile}.result"
  run_test "#{$bin}gt maxpairs -l 14 -ii sfxidx"
  run "cmp -s #{resultfile} #{$last_stdout}"
end

def checkmaxpairsiwithquery(reffile,queryfile)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)
  idxname=reffile + "-idx"
  run_test "#{$bin}gt suffixerator -algbds 3 40 120 -db " +
           "#{reffilepath} -indexname #{idxname} -dna -suf -tis -lcp -pl"
  run_test "#{$bin}gt maxpairs -l 15 -ii #{idxname} -q #{queryfilepath}"
  run "sort #{$last_stdout}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/mkvtree.x -indexname mkv-idx " +
  #    "-allout -v -pl -dna -db #{reffilepath}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/vmatch-mini.x 15 mkv-idx " +
  #    "#{queryfilepath}"
  #run "sed -e '/^#/d' #{$last_stdout} | sort"
  run "cmp -s #{$gttestdata}/maxpairs-result/#{reffile}-#{queryfile}.result #{$last_stdout}"
end

maxpairstestfiles.each do |reffile|
  Name "gt maxpairs #{reffile}"
  Keywords "gt_maxpairs gttestdata"
  Test do
    checkmaxpairs(reffile)
  end
  maxpairstestfiles.each do |queryfile|
    if reffile != queryfile
      Name "gt maxpairs #{reffile} versus #{queryfile}"
      Keywords "gt_maxpairs gttestdata"
      Test do
        checkmaxpairsiwithquery(reffile,queryfile)
      end
    end
  end
end
exit 1

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt greedyfwdmat #{reffile} #{queryfile}"
      Keywords "gt_greedyfwdmat small"
      Test do
        createandcheckgreedyfwdmat("#{$testdata}/#{reffile}",
                                   "#{$testdata}/#{queryfile}")
        checktagerator("#{$testdata}/#{reffile}",
                       "#{$testdata}/#{queryfile}")
        run "rm -f sfx.* fmi.* pck.*"
      end
    end
  end
end

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt idxlocali #{reffile} #{queryfile}"
      Keywords "gt_idxlocali"
      Test do
        run("#{$bin}gt packedindex mkindex -ssp -tis -indexname pck -db " +
            "#{$testdata}/#{reffile} -sprank -dna -pl -bsize 10 " +
            "-locfreq 32 -dir rev", 
            :maxtime => 100)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck " +
                 "-q #{$testdata}/#{queryfile}",
                 :maxtime => 100)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck -online " +
                 "-q #{$testdata}/#{queryfile}",
                 :maxtime => 100)
        run_test "#{$bin}gt suffixerator -indexname sfx -ssp -tis -suf -dna " +
                 "-v -db #{$testdata}/#{reffile}"
        run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx " +
                 "-q #{$testdata}/#{queryfile}", 
                 :maxtime => 100)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx -online " +
                 "-q #{$testdata}/#{queryfile}",
                 :maxtime => 100)
      end
    end
  end
end

allfiles.each do |reffile|
  Name "gt packedindex #{reffile}"
  Keywords "gt_packedindex small"
  Test do
    run_test("#{$bin}gt packedindex mkindex -tis -ssp -indexname pck " +
             "-sprank -db #{$testdata}/#{reffile} -dna -pl -bsize 10 " +
             " -locfreq 32 -dir rev", 
             :maxtime => 1200)
  end
end

Name "gt greedyfwdmat at1MB U8"
Keywords "gt_greedyfwdmat gttestdata"
Test do
  createandcheckgreedyfwdmat("#{$testdata}at1MB",
                             "#{$testdata}U89959_genomic.fas")
  run "rm -f sfx.* fmi.* pck.*"
end
