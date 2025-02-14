
#include "../include/utils.h"
#include <gtest/gtest.h>



TEST(TestCharBuffer, Constructor) {
    size_t buff_size { 100 };
    CharBuffer buff { buff_size };

    EXPECT_EQ(buff.buffer_size(), buff_size);
}


TEST(TestCharBuffer, Append) {
    size_t buff_size { 5 };
    CharBuffer buff { buff_size };

    buff.append('t');
    EXPECT_EQ(buff(0), 't');
    
    buff.append('h');
    buff.append('e');
    EXPECT_EQ(buff(0), 't');
    EXPECT_EQ(buff(1), 'h');
    EXPECT_EQ(buff(2), 'e');

    EXPECT_THROW(buff(6), std::out_of_range);
    EXPECT_THROW(buff(5), std::out_of_range);
    EXPECT_THROW(buff(-1), std::out_of_range);

    buff.append('r');
    buff.append('e');
    EXPECT_THROW(buff.append('e'), std::out_of_range); 
}


TEST(TestCharBuffer, Reset) {
    size_t buff_size { 5 };
    CharBuffer buff { buff_size };

    buff.append('t');
    EXPECT_EQ(buff(0), 't');

    buff.reset();
    EXPECT_THROW(buff(0), std::out_of_range);
}


TEST(TestCharBuffer, Data) {
    size_t buff_size { 5 };
    CharBuffer buff { buff_size };
    buff.append('t');
    buff.append('h');
    buff.append('e');

    EXPECT_EQ(buff.data()[0], 't');

    std::string s { buff.data() };
    EXPECT_EQ(s, "the");

    buff.append('e');
    s = buff.data();
    EXPECT_EQ(s, "thee");

}


TEST(TestStringRecord, ConstructorDelim) {
    CharBuffer buf { 10 };
    char s[] { "the\tcat\n" };
    StringRecord record { '\t' };

    record.update_str(s);

    EXPECT_TRUE(record.next_field());
    EXPECT_EQ(static_cast<std::string>(record.data()), "the");

    EXPECT_TRUE(record.next_field());
    EXPECT_EQ(static_cast<std::string>(record.data()), "cat");

    EXPECT_FALSE(record.next_field());

    EXPECT_EQ(record.size(), 8);

    record.reset();

    EXPECT_EQ(record.size(), 0);
}

TEST(TestStringRecord, SingleArgConstructor) {
    StringRecord record { '\t' };

    EXPECT_EQ(record.size(), 0);
    EXPECT_FALSE(record.next_field());

}
