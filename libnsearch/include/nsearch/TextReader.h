#pragma once

#include <string>
#include <iostream>
#include <fstream>

class TextReader {
public:
  virtual size_t NumBytesRead() const = 0;
  virtual size_t NumBytesTotal() const = 0;

  virtual bool EndOfFile() const = 0;

  virtual void operator>>( std::string &str ) = 0;
};

class TextStreamReader : public TextReader {
public:
  TextStreamReader( std::istream &is );

  size_t NumBytesRead() const;
  size_t NumBytesTotal() const;

  bool EndOfFile() const;

  void operator>>( std::string &str );

private:
  std::istream &mInput;
  std::streampos mTotalBytes;
};

/*
 * Reading huge fastq files needs to be fast
 */
#define TEXT_FILE_READER_BUFFER_SIZE 32 * 1024

class TextFileReader : public TextReader {
public:
  TextFileReader( const std::string &fileName );
  ~TextFileReader();

  size_t NumBytesRead() const;
  size_t NumBytesTotal() const;

  bool EndOfFile() const;

  void operator>>( std::string &str );

private:
  void NextBuffer();

  int mFd;

  size_t mBufferPos, mBufferSize;
  char mBuffer[ TEXT_FILE_READER_BUFFER_SIZE ];
  off_t mTotalBytes;
};
