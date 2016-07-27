#/usr/bin/ruby -w

unless ARGV.length == 5
  puts "This script modifies tab delimited file by replacing names in one column with new ones based on input file. Input arguments are (1) tab-delimited file with names to swap and replacement names (2) column number, starting with 0, containing term to search for, (3) column number containg term that will replace searched term (4) tab-delimited file with data fields to be converted (5) column in data file to replace"
  exit(1)
end

puts "Setting up replacement table to replace term in column #{ARGV[1]} with #{ARGV[2]} of input file #{ARGV[0]}"

list = Hash.new
File.open(ARGV[0]) do |file|
  file.each do |line|
    line.chomp!
    list[line.split("\t")[ARGV[1].to_i]] = line.split("\t")[ARGV[2].to_i]
  end
end

puts "Loading search list and adding target table term and writing to outfile"

outfile1 = File.new("#{ARGV[3]}"'.renamed',"w")
File.open(ARGV[3]) do |file|
  file.each do |line|
    line.chomp!
    i = ARGV[4].to_i
    temp = line.split("\t")
    temp[i] = list[temp[i]]
    outfile1.puts "#{temp.join("\t")}"
  end
end

outfile1.close

