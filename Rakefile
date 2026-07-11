# ----------------------------------------------------------------------------
# Rakefile per Alglin
# ----------------------------------------------------------------------------
# Utilizza il nuovo CMakeLists.txt moderno (FetchContent per le dipendenze).
# Non richiede sottomoduli git né cmake_utils esterni.
# ----------------------------------------------------------------------------

require 'colorize'
require 'fileutils'
require 'rake/clean'
require 'shellwords'

# Evita di cancellare file "core" (Eigen)
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

# ----------------------------------------------------------------------------
# Configurazione esterna (opzionale)
# ----------------------------------------------------------------------------
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../../Rakefile_configure.rb'
else
  # Configurazione di default
  COMPILE_DEBUG      = false        # true -> Debug, false -> Release
  COMPILE_DYNAMIC    = false        # true -> build anche shared library
  COMPILE_EXECUTABLE = true         # true -> compila tests
end

# ----------------------------------------------------------------------------
# Rilevamento OS
# ----------------------------------------------------------------------------
case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux|cygwin/
  OS = :linux
when /msys/
  OS = :mingw
else
  OS = :win
end

# Backend BLAS/LAPACK di default. Su macOS usa Accelerate, come LapackWrapper.
if OS == :mac
  DEFAULT_BACKEND = 'Accelerate'
else
  DEFAULT_BACKEND = 'OpenBLAS'
end

# Backend accettati dal Rakefile. Il valore canonico e' quello passato a
# LAPACK_WRAPPER_BACKEND, cioe' la stessa opzione usata dal Rakefile di
# LapackWrapper.
LAPACK_BACKEND_ALIASES = {
  'accelerate'       => 'Accelerate',
  'appleaccelerate'  => 'Accelerate',
  'apple_accelerate' => 'Accelerate',
  'openblas'         => 'OpenBLAS',
  'open_blas'        => 'OpenBLAS',
  'mkl'              => 'MKL',
  'blasfeo'          => 'BLASFEO'
}.freeze

# Il CMakeLists di Alglin mantiene proprie opzioni ALGLIN_USE_* per
# risolvere le librerie BLAS/LAPACK linkate direttamente da Alglin_Static e dai
# test. Teniamole sincronizzate con il backend scelto per LapackWrapper.
ALGLIN_BACKEND_FLAGS = {
  'Accelerate' => 'ALGLIN_USE_ACCELERATE',
  'OpenBLAS'   => 'ALGLIN_USE_OPENBLAS',
  'MKL'        => 'ALGLIN_USE_MKL'
}.freeze

def normalize_lapack_backend(backend)
  key = backend.to_s.strip.gsub(/[-\s]/, '_').downcase
  canonical = LAPACK_BACKEND_ALIASES[key]
  unless canonical
    valid = LAPACK_BACKEND_ALIASES.values.uniq.join(', ')
    raise ArgumentError, "Backend BLAS/LAPACK non supportato: #{backend.inspect}. Usa uno tra: #{valid}"
  end
  canonical
end

def lapack_backend_cmake_args(backend)
  backend = normalize_lapack_backend(backend)

  args = [
    "-DLAPACK_WRAPPER_BACKEND=#{backend}"
  ]

  # Passa solo l'opzione ALGLIN_USE_* effettivamente necessaria.  Non passare
  # opzioni legacy LAPACK_WRAPPER_USE_* né flag *_OFF: le revisioni moderne di
  # LapackWrapper usano LAPACK_WRAPPER_BACKEND, e CMake segnala come warning le
  # variabili di command line che non vengono lette.
  alglin_flag = ALGLIN_BACKEND_FLAGS[backend]
  args << "-D#{alglin_flag}:BOOL=ON" if alglin_flag

  args
end

# ----------------------------------------------------------------------------
# Helper: tipo di build
# ----------------------------------------------------------------------------
def build_type
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

# ----------------------------------------------------------------------------
# Helper: percorso radice del progetto
# ----------------------------------------------------------------------------
def project_root
  File.expand_path(__dir__)
end

# ----------------------------------------------------------------------------
# Helper: directory di build
# ----------------------------------------------------------------------------
def build_dir
  File.join(project_root, 'build')
end

# ----------------------------------------------------------------------------
# Helper: directory di installazione (lib/)
# ----------------------------------------------------------------------------
def install_prefix
  File.join(project_root, 'lib')
end

# ----------------------------------------------------------------------------
# Helper: cleanup del layout librerie locale
# ----------------------------------------------------------------------------
def cleanup_library_layout
  lib_root = File.join(project_root, 'lib')
  lib_dir  = File.join(lib_root, 'lib')

  FileUtils.mkdir_p lib_dir
  FileUtils.mkdir_p File.join(project_root, 'lib3rd', 'include')
  FileUtils.mkdir_p File.join(project_root, 'lib3rd', 'lib')

  # Residui generati da vecchie cache/install rules: il root lib/ non deve
  # contenere archivi, e non deve comparire lib/Users/...
  FileUtils.rm_f Dir.glob(File.join(lib_root, 'libAlglin*'))
  FileUtils.rm_rf File.join(lib_root, 'Users')

  # lib/lib deve contenere solo gli artifact di Alglin.  Le librerie delle
  # dipendenze non vanno installate qui.
  patterns = %w[
    libUtilsLite*
    liblapack_wrapper*
    libLapackWrapper*
    libHSL*
  ]
  patterns.each do |pattern|
    Dir.glob(File.join(lib_dir, pattern)).each { |path| FileUtils.rm_f(path) }
  end
end

# ----------------------------------------------------------------------------
# Generatore CMake (Ninja è il più veloce e cross‑platform)
# ----------------------------------------------------------------------------
def cmake_generator
  'Ninja'
end

# ----------------------------------------------------------------------------
# Comando CMake per la configurazione
# ----------------------------------------------------------------------------
def cmake_configure_command(backend, enable_tests = false)
  backend = normalize_lapack_backend(backend)
  enable_tests_value = enable_tests ? 'ON' : 'OFF'

  args = [
    "-G", cmake_generator,
    "-B", build_dir,
    "-DCMAKE_BUILD_TYPE=#{build_type}",
    "-DCMAKE_INSTALL_PREFIX=#{install_prefix}",
    "-DCMAKE_INSTALL_LIBDIR=lib",
    "-DCMAKE_INSTALL_INCLUDEDIR=include",
    "-DUTILS_BUILD_SHARED=#{COMPILE_DYNAMIC ? 'ON' : 'OFF'}",
    "-DUTILS_ENABLE_TESTS=#{enable_tests_value}",
    "-DALGLIN_ENABLE_TESTS=#{enable_tests_value}",
    "-DUTILS_INSTALL_LOCAL=OFF"  # installa in lib/ con layout lib/include + lib/lib
  ]

  args.concat lapack_backend_cmake_args(backend)

  # Se si usa OpenBLAS su Windows, aggiungi l'opzione BITS prima del path sorgente.
  args << "-DBITS=64" if backend == 'OpenBLAS' && OS == :win

  args << project_root
  args.shelljoin
end

# ----------------------------------------------------------------------------
# Comando CMake per la build
# ----------------------------------------------------------------------------
def cmake_build_command(target = nil)
  cmd = ["--build", build_dir, "--config", build_type]
  cmd += ["--target", target] if target
  cmd << "--parallel"
  cmd.join(' ')
end

# ----------------------------------------------------------------------------
# Modifica del file di configurazione per disabilitare OpenMP
# (mantenuto per compatibilità, ma il nuovo CMake può gestirlo via opzioni)
# ----------------------------------------------------------------------------
def change_config_file
  file_path = File.join(project_root, 'src', 'Alglin_Config.hh')
  tmpl_path = file_path + '.tmpl'
  if File.exist?(tmpl_path)
    text = File.read(tmpl_path)
    text.gsub!('@@ALGIN_OPENMP@@', '#define ALGLIN_DO_NOT_USE_OPENMP 1')
    File.open(file_path, 'w+') { |f| f << text }
    puts "Aggiornato Alglin_Config.hh".yellow
  else
    puts "Nessun template trovato per Alglin_Config.hh".yellow
  end
end

# ----------------------------------------------------------------------------
# Task
# ----------------------------------------------------------------------------
desc "default task --> build"
task :default => :build

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end

desc "Configure CMake con un backend BLAS/LAPACK (Accelerate, OpenBLAS, MKL)"
task :configure, [:backend] do |t, args|
  backend = normalize_lapack_backend(args[:backend] || DEFAULT_BACKEND)
  puts "Configurazione CMake in #{build_dir} con backend #{backend}...".green
  change_config_file
  FileUtils.rm_rf build_dir if File.exist?(build_dir)
  sh "cmake " + cmake_configure_command(backend, false)
end

desc "Configure CMake con i test Alglin abilitati"
task :configure_tests, [:backend] do |t, args|
  backend = normalize_lapack_backend(args[:backend] || DEFAULT_BACKEND)
  puts "Configurazione CMake in #{build_dir} con backend #{backend} e test abilitati...".green
  change_config_file
  FileUtils.rm_rf build_dir if File.exist?(build_dir)
  sh "cmake " + cmake_configure_command(backend, true)
end

desc "Compila e installa la libreria (backend di default)"
task :build, [:backend] => :configure do |t, args|
  puts "Compilazione e installazione...".green
  sh "cmake " + cmake_build_command('install')
  cleanup_library_layout
end

desc "Compila tutti i target (inclusi i test, se abilitati dalla configurazione)"
task :compile_all, [:backend] => :configure do |t, args|
  puts "Compilazione di tutti i target...".green
  sh "cmake " + cmake_build_command
end

desc "Compila i test Alglin"
task :build_tests, [:backend] => :configure_tests do |t, args|
  puts "Compilazione dei test Alglin...".green
  sh "cmake " + cmake_build_command('Alglin_all_tests')
end

desc "Esegue i test Alglin (configura e compila prima i test)"
task :test, [:backend] => :build_tests do |t, args|
  puts "Esecuzione dei test Alglin...".green
  Dir.chdir(build_dir) do
    sh "ctest -C #{build_type} --output-on-failure"
  end
  cleanup_library_layout
end

desc "Configura, compila ed esegue i test Alglin"
task :run, [:backend] => :test do
  # stesso di test
end

desc "Pulisce la directory di build (mantiene lib/ e lib3rd/)"
task :clean_build do
  puts "Rimozione directory build...".green
  FileUtils.rm_rf build_dir
end

desc "Pulizia completa (rimuove build, lib, lib3rd, bin, object files)"
task :clean do
  puts "Pulizia completa...".green
  FileUtils.rm_rf build_dir
  FileUtils.rm_rf File.join(project_root, 'lib')
  FileUtils.rm_rf File.join(project_root, 'lib3rd')
  FileUtils.rm_rf File.join(project_root, 'bin')
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.o'))
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.obj'))
end

desc "Crea pacchetti con CPack"
task :cpack do
  puts "Creazione pacchetti...".green
  Dir.chdir(build_dir) do
    sh "cpack -C #{build_type} CPackConfig.cmake"
    sh "cpack -C #{build_type} CPackSourceConfig.cmake"
  end
end

# ----------------------------------------------------------------------------
# Task di compatibilità (come nel Rakefile originale)
# ----------------------------------------------------------------------------
task build_osx:   :build
task build_linux: :build
task build_mingw: :build
task build_win:   :build

task clean_osx:   :clean
task clean_linux: :clean
task clean_mingw: :clean
task clean_win:   :clean

# ----------------------------------------------------------------------------
# Pulizia automatica con Rake::Clean
# ----------------------------------------------------------------------------
CLOBBER.include(
  build_dir,
  File.join(project_root, 'lib'),
  File.join(project_root, 'lib3rd'),
  File.join(project_root, 'bin')
)
