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

cmakeversion = %x( cmake --version ).scan(/\d+\.\d+/).last
if cmakeversion >= "3.12" then
  PARALLEL = '--parallel 8 '
else
  PARALLEL = ''
end

TESTS = [
  "Simplex-Test1",
  "Simplex-Test2",
  "Simplex-Test3",
  "Simplex-Test4",
  "test0-FD",
  "test1-small-factorization",
  "test2-Threads",
  "test3-Timing",
  "test4-KKT",
  "test5-ABD-Diaz",
  "test6-ABD-Block",
  "test7-BorderedCR",
  "test8-Cinterface",
  "test12-BandedMatrix",
  "test13-BFGS",
  "test14-BLOCKTRID",
  "test15-EIGS"
]

desc "run tests"
task :run do
  TESTS.each do |cmd|
    sh "./bin/#{cmd}"
  end
end

task :run_win do
  TESTS.each do |cmd|
    sh "bin\\Release\\#{cmd}.exe"
  end
end

desc "build lib"
task :build2 do
  sh "make; make --jobs=8 all"
end

desc "build lib"
task :build, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
  FileUtils.rm_rf 'build'
  FileUtils.mkdir 'build'
  FileUtils.cd    'build'
  sh 'cmake -D'+args.lapack+':=true ..'
  sh 'make --jobs=8 install'
  FileUtils.cd '..'
end

desc "build lib"
task :build_exe, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
  FileUtils.cd 'build'
  sh 'cmake -D'+args.lapack+':=true -DBUILD_EXECUTABLE:=true ..'
  sh 'make --jobs=8 install'
  FileUtils.cd '..'
end

task :build_osx, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
  Rake::Task[:osx_3rd].invoke(args.lapack)
  Rake::Task[:build].invoke(args.lapack)
end

task :build_linux, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
  Rake::Task[:linux_3rd].invoke(args.lapack)
  Rake::Task[:build].invoke(args.lapack)
end

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text= File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:year, :bits] => [:win_3rd] do |t, args|
  args.with_defaults( :year => "2017", :bits => "x64" )

  puts "\n\nBUILD\n\n".green

  cmd = "set path=%path%;lib3rd\\lib;lib3rd\\dll;"

  FileUtils.rm_f 'src/AlglinSuperLU.hh'
  FileUtils.cp   'src/AlglinSuperLU.hh.tmpl', 'src/AlglinSuperLU.hh'

  ChangeOnFile(
    'src/AlglinSuperLU.hh',
    '@@VSYEARANDBITS@@',
    "vs#{args.year}_#{args.bits}"
  )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  # do not build executable
  #tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " + ' -DBUILD_EXECUTABLE=1 -DCMAKE_INSTALL_PREFIX:PATH=..\lib ..'
  tmp = " -DBITS=#{args.bits} -DYEAR=#{args.year} " +
        ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib ' +
        ' -DBUILD_EXECUTABLE:VAR=true ..'

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

  sh 'cmake --build . --config Release  --target ALL_BUILD '+PARALLEL
  FileUtils.mkdir_p "../lib"
  FileUtils.cp 'Release/Alglin.lib', "../lib/Alglin_vs#{args.year}_#{args.bits}.lib"
  sh 'cmake --build . --config Debug --target ALL_BUILD '+PARALLEL
  FileUtils.cp 'Debug/Alglin.lib', "../lib/Alglin_vs#{args.year}_#{args.bits}_debug.lib"

  FileUtils.cd '..'

end

##### desc 'install third parties for osx [lapack=LAPACK_WRAPPER_USE_ACCELERATE]'
##### task :osx_3rd, [:lapack] do |t, args|
#####   args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
#####   FileUtils.cd 'third_parties'
#####   sh "rake install_osx"
#####   FileUtils.cd '..'
#####   FileUtils.cd 'submodules'
#####   sh "rake build_osx[#{args.lapack}]"
#####   FileUtils.cd '..'
##### end
##### 
##### task :copy_3rd do
#####   FileUtils.mkdir_p "lib"
#####   FileUtils.mkdir_p "lib/lib"
#####   FileUtils.mkdir_p "lib/bin"
#####   FileUtils.mkdir_p "lib/dll"
#####   ['./submodules/LapackWrapper/lib',
#####    './submodules/LapackWrapper/lib3rd'].each do |path|
#####     FileUtils.cp_r "#{path}/lib",     'lib3rd' if Dir.exist?("#{path}/lib")
#####     FileUtils.cp_r "#{path}/dll",     'lib3rd' if Dir.exist?("#{path}/dll")
#####     FileUtils.cp_r "#{path}/bin",     'lib3rd' if Dir.exist?("#{path}/bin")
#####     FileUtils.cp_r "#{path}/include", 'lib3rd' if Dir.exist?("#{path}/include")
#####   end
#####   FileUtils.cd '..'
##### end

desc 'install third parties for osx [lapack=LAPACK_WRAPPER_USE_ACCELERATE]'
task :osx_3rd, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES\n\n".green
  sh "rake build_osx[#{args.lapack}]"
  FileUtils.cd '../third_parties'
  puts "\n\nTHIRD PARTIES\n\n".green
  sh "rake install_osx[#{args.lapack}]"
  FileUtils.cd '..'
end

desc 'install third parties for linux [lapack=LAPACK_WRAPPER_USE_OPENBLAS]'
task :linux_3rd, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES\n\n".green
  sh "rake build_linux[#{args.lapack}]"
  FileUtils.cd '../third_parties'
  puts "\n\nTHIRD PARTIES\n\n".green
  sh "rake install_linux[#{args.lapack}]"
  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :win_3rd, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
  )
  FileUtils.cd 'submodules'
  puts "\n\nSUBMODULES\n\n".green
  sh "rake build_win[#{args.year},#{args.bits},#{args.lapack}]"
  FileUtils.cd '../third_parties'
  puts "\n\nTHIRD PARTIES\n\n".green
  sh "rake install_win[#{args.year},#{args.bits},#{args.lapack}]"
  FileUtils.cd '..'
end

desc "clean for osx"
task :clean_osx do
  sh "make clean"
  FileUtils.cd 'third_parties'
  sh "rake clean_osx"
  FileUtils.cd '../submodules'
  sh "rake clean_osx"
  FileUtils.cd '..'
end

desc "clean for linux"
task :clean_linux do
  sh "make clean"
  FileUtils.cd 'third_parties'
  sh "rake clean_linux"
  FileUtils.cd '../submodules'
  sh "rake clean_linux"
  FileUtils.cd '..'
end

desc "clean for windows"
task :clean_win do
  FileUtils.rm_rf 'vs_*'
  FileUtils.cd 'third_parties'
  sh "rake clean_win"
  FileUtils.cd '../submodules'
  sh "rake clean_win"
  FileUtils.cd '..'
end
