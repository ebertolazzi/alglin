require_relative "./cmake_utils/Rakefile_common.rb"

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text = File.read file+".tmpl"
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

#task :mkl, [:year, :bits] do |t, args|
#  args.with_defaults(:year => "2017", :bits => "x64" )
#  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
#end

desc "compile for OSX [default lapack=LAPACK_WRAPPER_USE_ACCELERATE]"
task :build_osx, [:lapack] do |t, args|
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
  Rake::Task[:osx_3rd].invoke(args.lapack)
  Rake::Task[:build_common].invoke(args.lapack)
end

desc "compile for LINUX [default lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_linux, [:lapack] do |t, args|
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
  Rake::Task[:linux_3rd].invoke(args.lapack)
  Rake::Task[:build_common].invoke(args.lapack)
end

desc "compile for MINGW [default lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_mingw, [:lapack] do |t, args|
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
  Rake::Task[:mingw_3rd].invoke(args.lapack)
  Rake::Task[:build_common].invoke(args.lapack)
end

desc "compile for OSX [default lapack=LAPACK_WRAPPER_USE_ACCELERATE]"
task :build_common, [:lapack] do |t, args|
  ChangeOnFile(
    'src/Alglin_Config.hh',
    '@@ALGIN_OPENMP@@',
    "#define ALGLIN_DO_NOT_USE_OPENMP 1"
  )

  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  puts "run CMAKE for ALGLIN".yellow
  sh "cmake -G Ninja " + cmd_cmake_build() + ' -DUTILS_USE_OPENMP=OFF -D' + args.lapack + ':VAR=ON ..'

  puts "compile with CMAKE for ALGLIN".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017 bits=x64]"
task :build_win, [:lapack] do |t, args|

  # check architecture
  case `where cl.exe`.chop
  when /x64\\cl\.exe/
    VS_ARCH = 'x64'
  when /amd64\\cl\.exe/
    VS_ARCH = 'x64'
  when /bin\\cl\.exe/
    VS_ARCH = 'x86'
  else
    raise RuntimeError, "Cannot determine architecture for Visual Studio".red
  end

  ChangeOnFile(
    'src/Alglin_Config.hh',
    '@@ALGIN_OPENMP@@',
    "#define ALGLIN_DO_NOT_USE_OPENMP 1"
  )

  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'

  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )

  Rake::Task[:win_3rd].invoke(args.lapack)

  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+VS_ARCH
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for ALGLIN".yellow
  sh "cmake -G Ninja -DBITS:VAR=#{VS_ARCH} " + cmd_cmake_build() + ' -DUTILS_USE_OPENMP=OFF -D' + args.lapack + ':VAR=ON ..'

  puts "compile with CMAKE for ALGLIN".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

desc 'install third parties for osx [lapack=LAPACK_WRAPPER_USE_ACCELERATE]'
task :osx_3rd, [:lapack] do |t, args|
  if File.directory?('../LapackWrapper/lib') then
    FileUtils.cp_r '../LapackWrapper/lib/.', 'lib3rd'
    FileUtils.cp_r '../LapackWrapper/lib3rd/.', 'lib3rd' if File.directory?('../LapackWrapper/lib3rd')
  else
    args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )
    FileUtils.cd 'submodules'
      puts "\n\nSUBMODULES (for ALGLIN)\n\n".green
      sh "rake build_osx[#{args.lapack}]"
    FileUtils.cd '..'
  end
end

desc 'install third parties for linux [lapack=LAPACK_WRAPPER_USE_OPENBLAS]'
task :linux_3rd, [:lapack] do |t, args|
  if File.directory?('../LapackWrapper/lib') then
    FileUtils.cp_r '../LapackWrapper/lib/.', 'lib3rd'
    FileUtils.cp_r '../LapackWrapper/lib3rd/.', 'lib3rd' if File.directory?('../LapackWrapper/lib3rd')
  else
    args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
    FileUtils.cd 'submodules'
      puts "\n\nSUBMODULES (for ALGLIN)\n\n".green
      sh "rake build_linux[#{args.lapack}]"
    FileUtils.cd '..'
  end
end

desc 'install third parties for MINGW [lapack=LAPACK_WRAPPER_USE_OPENBLAS]'
task :mingw_3rd, [:lapack] do |t, args|
  if File.directory?('../LapackWrapper/lib') then
    FileUtils.cp_r '../LapackWrapper/lib/.', 'lib3rd'
    FileUtils.cp_r '../LapackWrapper/lib3rd/.', 'lib3rd' if File.directory?('../LapackWrapper/lib3rd')
  else
    args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
    FileUtils.cd 'submodules'
      puts "\n\nSUBMODULES (for ALGLIN)\n\n".green
      sh "rake build_mingw[#{args.lapack}]"
    FileUtils.cd '..'
  end
end

desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :win_3rd, [:lapack] do |t, args|
  if File.directory?('../LapackWrapper/lib') then
    FileUtils.cp_r '../LapackWrapper/lib/.', 'lib3rd'
    FileUtils.cp_r '../LapackWrapper/lib3rd/.', 'lib3rd' if File.directory?('../LapackWrapper/lib3rd')
  else
    args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_OPENBLAS" )
    FileUtils.cd 'submodules'
      puts "\n\nSUBMODULES (for ALGLIN)\n\n".green
      sh "rake build_win[#{args.lapack}]"
    FileUtils.cd '..'
  end
end

desc "clean common"
task :clean_common do
  FileUtils.rm_rf 'bin'
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  #sh "make clean"
  FileUtils.cd 'submodules'
  sh "rake clean"
  FileUtils.cd '..'
end

task :clean_osx   => :clean_common do end
task :clean_linux => :clean_common do end
task :clean_mingw => :clean_common do end
task :clean_win   => :clean_common do end

task :cppcheck do
  FileUtils.rm_rf   'lib'
  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'
  sh 'cmake -DCMAKE_EXPORT_COMPILE_COMMAND=ON ..'
  sh 'cppcheck --project=compile_commands.json'
end

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for ROOTS".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
