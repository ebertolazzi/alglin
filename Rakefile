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
task :run do
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

task :run_win, [:year, :bits] do |t, args|
  args.with_defaults( :year => "2017", :bits => "x64" )
  sh "runtest.bat"
end

desc "build lib"
task :build  do
  sh "make"
end

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text= File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:year, :bits, :lapack, :thread] do |t, args|
  args.with_defaults( :year   => "2017",
                      :bits   => "x64",
                      :lapack => "ALGLIN_USE_LAPACK",
                      #:lapack => "ALGLIN_USE_LAPACK2",
                      #:lapack => "ALGLIN_USE_OPENBLAS",
                      :thread => "ALGLIN_USE_THREAD" )

  cmd = "set path=%path%;lib3rd\\lib;lib3rd\\dll;"

  FileUtils.rm_f 'src/AlglinConfig.hh'
  FileUtils.cp   'src/AlglinConfig.hh.tmpl', 'src/AlglinConfig.hh'
  FileUtils.rm_f 'src/AlglinSuperLU.hh'
  FileUtils.cp   'src/AlglinSuperLU.hh.tmpl', 'src/AlglinSuperLU.hh'

  ChangeOnFile( 'src/AlglinConfig.hh',
                '@@ALGLIN_USE@@',
                "#define #{args.lapack} 1" )
  ChangeOnFile( 'src/AlglinConfig.hh',
                '@@ALGLIN_THREAD@@',
                "#define #{args.thread} 1" )
  ChangeOnFile( 'src/AlglinConfig.hh',
                '@@ALGLIN_NOSYSTEM_OPENBLAS@@',
                "#define ALGLIN_DO_NOT_USE_SYSTEM_OPENBLAS 1" )
  ChangeOnFile( 'src/AlglinSuperLU.hh',
                '@@VSYEARANDBITS@@',
                "vs#{args.year}_#{args.bits}" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  # do not build executable
  #tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " + ' -DBUILD_EXECUTABLE=1 -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..'
  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " + ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..'

  win32_64 = ''
  case args.bits
  when /x64/
    win32_64 = ' Win64'
  end

  case args.year
  when "2010"
    sh 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
  when "2012"
    sh 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
  when "2013"
    sh 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
  when "2015"
    sh 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
  when "2017"
    sh 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
  else
    puts "Visual Studio year #{year} not supported!\n";
  end

  sh 'cmake --build . --config Release  --target ALL_BUILD'
  FileUtils.mkdir_p "../lib"
  FileUtils.cp 'Release/Alglin.lib', "../lib/Alglin_vs#{args.year}_#{args.bits}.lib"  
  sh 'cmake --build . --config Debug --target ALL_BUILD'
  FileUtils.cp 'Debug/Alglin.lib', "../lib/Alglin_vs#{args.year}_#{args.bits}_debug.lib"

  FileUtils.cd '..'

end
