add_targets('mySong.txt')

rule sing_a_song:
  input: EF("[PP:songs][P:song]")
  output: touch(T("mySong.txt"))
  params: song=P('song'), name=P('name'), pp=PP('songs')
  run:
    correct = {
      "P/Peter/mySong.txt": ["Peter", "input", "songs", "1.txt"],
      "P/Paul/mySong.txt": ["Paul", "input", "songs", "2.txt"],
      "P/Mary/mySong.txt": ["Mary", "input", "songs", "3.txt"]
    }	
    f = input[0].split("/")
    assert output[0] in correct
    assert params.name == correct[output[0]][0]
    assert params.song == correct[output[0]][3]
    assert f[-1] == correct[output[0]][3]
    assert f[-2] == correct[output[0]][2]
    assert f[-3] == correct[output[0]][1]

         
