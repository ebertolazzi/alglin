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

require_relative "./Rakefile_common.rb"

task :default => [:build]

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
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

desc "compile for OSX [default lapack=LAPACK_WRAPPER_USE_ACCELERATE]"
task :build_osx, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )

  Rake::Task[:osx_3rd].invoke(args.lapack)

  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  cmake_cmd = 'cmake -DBUILD_SHARED:VAR='
  if COMPILE_DYNAMIC then
    cmake_cmd += 'true '
  else
    cmake_cmd += 'false '
  end
  cmake_cmd += '-D' + args.lapack + ':VAR=true '
  if COMPILE_EXECUTABLE then
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=false '
  end

  if COMPILE_DEBUG then
    sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end
  sh cmake_cmd + '-DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --config Release --target install '+PARALLEL
  FileUtils.cd '..'

end

desc "compile for LINUX [default lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_linux, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )

  Rake::Task[:linux_3rd].invoke(args.lapack)

  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  cmake_cmd = 'cmake -DBUILD_SHARED:VAR='
  if COMPILE_DYNAMIC then
    cmake_cmd += 'true '
  else
    cmake_cmd += 'false '
  end
  cmake_cmd += '-D' + args.lapack + ':VAR=true '
  if COMPILE_EXECUTABLE then
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += '-DBUILD_EXECUTABLE:VAR=false '
  end

  if COMPILE_DEBUG then
    sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end
  sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --config Release --target install '+PARALLEL
  FileUtils.cd '..'

end

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text= File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
  )

  Rake::Task[:win_3rd].invoke(args.year,args.bits,args.lapack)

  cmd = "set path=%path%;lib3rd\\lib;lib3rd\\dll;"
  dir = "vs_#{args.year}_#{args.bits}"

  # do not build executable
  cmake_cmd = win_vs(args.bits,args.year)
  cmake_cmd += " -DBITS=#{args.bits} -DYEAR=#{args.year} " +
               ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib '

  if COMPILE_EXECUTABLE then
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=true '
  else
    cmake_cmd += ' -DBUILD_EXECUTABLE:VAR=false '
  end
  if COMPILE_DYNAMIC then
    cmake_cmd += ' -DBUILD_SHARED:VAR=true '
  else
    cmake_cmd += ' -DBUILD_SHARED:VAR=false '
  end

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  if COMPILE_DEBUG then
    FileUtils.rm_rf   dir
    FileUtils.mkdir_p dir
    FileUtils.cd      dir

    sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --clean-first --config Debug --target install '+PARALLEL
    FileUtils.cd '..'
  end

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --clean-first --config Release  --target install '+PARALLEL
  FileUtils.cd '..'

end

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
