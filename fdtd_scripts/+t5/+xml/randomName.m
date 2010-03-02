function filename = randomName(prefix, suffix, numRandomChars)

filename = [prefix, char('a'+randint(1,numRandomChars,[0,25])), suffix];
