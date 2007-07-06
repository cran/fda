
Dir["*/*.[rR]d"].each do |path|
  puts path
  file = File.new(path, "r")

  update = file.read.gsub(/\\seealso\{(.*?)\}/smi) do
    "\\seealso{\n" +
    $1.gsub("\n","").split(/[\, ]+/).map{|f| "\\code{\\link{#{f}}}"}.join(", \n") +
    "\n}"
  end
  
  file.close

  file = File.new(path, "w")
  file.write update
  file.close
  
end
