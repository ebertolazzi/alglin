#
#
#

%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

task :default => [:build]

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
end

desc "run tests"
task :run  do
  sh "./bin/test0-FD"
  sh "./bin/test1-small-factorization"
  sh "./bin/test2-Threads"
  sh "./bin/test3-Timing"
  sh "./bin/test4-KKT"
  sh "./bin/test5-ABD-Diaz"
  sh "./bin/test6-ABD-Block"
  sh "./bin/test7-BorderedCR"
  sh "./bin/test8-Cinterface"
  sh "./bin/test9-Cinterface"
  sh "./bin/test12-BandedMatrix"
end

desc "build lib"
task :build  do
  sh "make"
end

task :build_win, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf dir
  FileUtils.mkdir dir
  FileUtils.cd    dir
  case args.year
  when "2010"
    sh 'cmake -G "Visual Studio 10 2010" -DCMAKE_INSTALL_PREFIX:PATH=..\libs ..'
  when "2012"
    sh 'cmake -G "Visual Studio 11 2012" -DCMAKE_INSTALL_PREFIX:PATH=..\libs ..'
  when "2013"
    sh 'cmake -G "Visual Studio 12 2013" -DCMAKE_INSTALL_PREFIX:PATH=..\libs ..'
  when "2015"
    sh 'cmake -G "Visual Studio 14 2015" -DCMAKE_INSTALL_PREFIX:PATH=..\libs ..'
  when "2017"
    sh 'cmake -G "Visual Studio 15 2017" -DCMAKE_INSTALL_PREFIX:PATH=..\libs ..'
  else
    puts "Visual Studio year #{args.year} not supported!\n";
  end

  sh 'cmake --build . --config Release  --target all'
  sh 'cmake --build . --config Debug --target all'

  FileUtils.cd '..'

end
