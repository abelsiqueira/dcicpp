#!/usr/bin/env ruby

special_variables = {
  "dciTrue"=>1,
  "dciFalse"=>0,
  "dciInf"=>1e20
}

puts File.open("src/parameters.cpp", "r").each_line.map { |line|
  if !(m = line.match(/paramMap\["(.*)"\]/)).nil?
    m[1]
  end
}.compact.map { |var|
  File.open("src/parameters.cpp", "r").each_line.map { |line|
    if !(m = line.match(/\b#{var} = (.*);/)).nil?
      if special_variables.include?(m[1])
        value = special_variables[m[1]]
      else
        value = m[1]
      end
      break "#{var} #{value}"
    end
  }
}.sort_by { |w| w.downcase }
