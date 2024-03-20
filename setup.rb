puts "Setup submodules"
system('git submodule foreach --recursive git submodule init')
system('git submodule foreach --recursive git submodule update')
system('git submodule foreach --recursive git submodule sync')
system('cd submodule/LapackWrapper; git checkout develop')
system('cd submodules/LapackWrapper/submodules/Utils; git checkout main')
system('cd submodules/LapackWrapper/submodules/Utils/cmake_utils; git checkout main')
system('cd submodules/cmake_utils; git checkout main')

puts ARGV

if ARGV.size() > 0 && ARGV[0] == "--last" then
  puts "\nUpdate submodules to last version"
  system('git submodule foreach --recursive git pull')
end
