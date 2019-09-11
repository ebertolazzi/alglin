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

require "./Rakefile_common.rb"

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

task :build_osx, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )

  Rake::Task[:osx_3rd].invoke(args.lapack)

  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  if COMPILE_DEBUG then
    sh 'cmake -D' + args.lapack + ':VAR=true -DBUILD_EXECUTABLE:VAR=true -DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end
  sh 'cmake -D' + args.lapack + ':VAR=true -DBUILD_EXECUTABLE:VAR=true -DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --config Release --target install '+PARALLEL
  FileUtils.cd '..'

end

task :build_linux, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )

  Rake::Task[:linux_3rd].invoke(args.lapack)

  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  if COMPILE_DEBUG then
    sh 'cmake -D' + args.lapack + ':VAR=true -DBUILD_EXECUTABLE:VAR=true -DCMAKE_BUILD_TYPE:VAR=Debug ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  end
  sh 'cmake -D' + args.lapack + ':VAR=true -DBUILD_EXECUTABLE:VAR=true -DCMAKE_BUILD_TYPE:VAR=Release ..'
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
               ' -DCMAKE_INSTALL_PREFIX:PATH=..\lib ' +
               ' -DBUILD_EXECUTABLE:VAR=true ..'

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
    sh 'cmake --build . --config Debug --target install '+PARALLEL
    FileUtils.cp_r './lib/dll', '../lib/' if Dir.exist?('./lib/dll')
    Dir['./lib/bin/*'].each do |f|
      FileUtils.cp f, '../lib/bin/'+args.bits+'/'+File.basename(f)
    end
    Dir['./lib/lib/*'].each do |f|
      if /\_static.*\.lib$/.match(f) then
        FileUtils.cp f, '../lib/lib/'+File.basename(f)
      else
        FileUtils.cp f, '../lib/dll/'+File.basename(f)
      end
    end
    FileUtils.cd '..'
  end

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  sh cmake_cmd + ' -DCMAKE_BUILD_TYPE:VAR=Release ..'
  sh 'cmake --build . --clean-first --config Release  --target install '+PARALLEL
  FileUtils.cp_r './lib/dll', '../lib/' if Dir.exist?('./lib/dll')
  Dir['./lib/bin/*'].each do |f|
    FileUtils.cp f, '../lib/bin/'+args.bits+'/'+File.basename(f)
  end
  Dir['./lib/lib/*'].each do |f|
    if /\_static.*\.lib$/.match(f) then
      FileUtils.cp f, '../lib/lib/'+File.basename(f)
    else
      FileUtils.cp f, '../lib/dll/'+File.basename(f)
    end
  end
  FileUtils.cp_r './lib/include', '../lib/' if Dir.exist?('./lib/include')
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
